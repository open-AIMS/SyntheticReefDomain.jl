# Assemble a synthetic domain in memory, and write it out as an ADRIA data package.

"""
In-memory synthetic domain (everything that gets written, plus the regional DHW signal and
cyclone intensities used for previews).

Fields:
- `config` : the validated `DomainConfig`
- `reefs` : one row per reef (see [`build_reefs`](@ref))
- `connectivity` : source x sink transfer probabilities (`n x n`)
- `cover` : initial cover, `5 x n`, fraction of habitable area
- `dhw` : climate label => DHW (`T x n x S_dhw`)
- `dhw_regional` : climate label => regional DHW signal (`T x S_dhw`)
- `waves` : wave Ub (`T x n x S_dhw`)
- `cyclone_mortality` : `T x n x 5 x S_cyc`
- `cyclone_intensity` : `T x S_cyc`, 0 in years without a cyclone
"""
struct SyntheticDomain
    config::DomainConfig
    reefs::DataFrame
    connectivity::Matrix{Float64}
    cover::Matrix{Float64}
    dhw::Dict{Int,Array{Float64,3}}
    dhw_regional::Dict{Int,Matrix{Float64}}
    waves::Array{Float64,3}
    cyclone_mortality::Array{Float64,4}
    cyclone_intensity::Matrix{Float64}
end

years(dom::SyntheticDomain) = years(dom.config)
climate_labels(dom::SyntheticDomain) = climate_labels(dom.config)

function Base.show(io::IO, ::MIME"text/plain", dom::SyntheticDomain)
    cfg = dom.config
    T = length(years(dom))
    n_events = count(>(0), dom.cyclone_intensity)
    println(io, "SyntheticDomain \"$(cfg.name)\"")
    println(io, "  reefs        : $(cfg.reefs.n) ($(cfg.reefs.layout) layout, " *
                "$(cfg.reefs.n_clusters) cluster(s))")
    println(io, "  timeframe    : $(cfg.timeframe[1])-$(cfg.timeframe[2]) ($T years)")
    println(io, "  connectivity : $(cfg.connectivity.topology), " *
                "$(count(>(0), dom.connectivity) - cfg.reefs.n) off-diagonal links")
    println(io, "  initial cover: $(round(100minimum(sum(dom.cover; dims=1)); digits=1))-" *
                "$(round(100maximum(sum(dom.cover; dims=1)); digits=1))% of habitable area")
    for label in climate_labels(dom)
        D = dom.dhw[label]
        first_dec = mean(@view D[1:min(10, T), :, :])
        last_dec = mean(@view D[max(1, T - 9):T, :, :])
        println(io, "  DHW (label $label): first decade $(round(first_dec; digits=1)) -> " *
                    "last decade $(round(last_dec; digits=1)) DHW, " *
                    "$(cfg.n_dhw_scenarios) scenarios")
    end
    print(io, "  cyclones     : ", cfg.cyclones.enabled ?
        "$(round(n_events / cfg.n_cyclone_scenarios; digits=1)) events per scenario on average, " *
        "$(cfg.n_cyclone_scenarios) scenarios" : "disabled")
    return nothing
end

"""
    build_domain(cfg::DomainConfig)::SyntheticDomain
    build_domain(path::AbstractString)::SyntheticDomain

Validate a configuration (or load it from a TOML file) and build the domain in memory,
without writing anything.
"""
build_domain(path::AbstractString)::SyntheticDomain = build_domain(load_config(path))
function build_domain(cfg::DomainConfig)::SyntheticDomain
    validate(cfg)
    reefs = build_reefs(cfg)
    connectivity = build_connectivity(cfg, reefs)
    cover = build_cover(cfg, reefs)
    dhw, regional = build_dhw(cfg, reefs)
    waves = build_waves(cfg, reefs)
    mortality, intensity = build_cyclones(cfg, nrow(reefs))
    return SyntheticDomain(
        cfg, reefs, connectivity, cover, dhw, regional, waves, mortality, intensity
    )
end

"""
    write_domain(dom::SyntheticDomain, outdir::AbstractString=joinpath("domains", dom.config.name); force=false)::String

Write `dom` as an ADRIA data package into `outdir` and return the path. ADRIA names the
domain after the folder, so the geopackage is called `<basename(outdir)>.gpkg`.

An existing non-empty folder is only replaced with `force=true`, and only if it already
holds a data package (`datapackage.json`), so an unrelated folder is never deleted.
"""
function write_domain(
    dom::SyntheticDomain, outdir::AbstractString=joinpath("domains", dom.config.name);
    force::Bool=false
)::String
    outdir = normpath(abspath(outdir))
    name = basename(outdir)
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", name) || throw(ArgumentError(
        "Output folder name '$name' must use only letters, digits, '_', '-' or '.'."
    ))
    name == dom.config.name || @info "Writing domain \"$(dom.config.name)\" as \"$name\" " *
                                     "(ADRIA names domains after their folder)."

    if isdir(outdir) && !isempty(readdir(outdir))
        force || throw(ArgumentError(
            "$outdir already exists and is not empty. Pass force=true to replace it."
        ))
        isfile(joinpath(outdir, "datapackage.json")) || throw(ArgumentError(
            "Refusing to replace $outdir: it does not contain a datapackage.json, so it " *
            "does not look like a domain written by this package."
        ))
        rm(outdir; recursive=true)
    end

    for sub in ("spatial", "connectivity", "DHWs", "waves", "cyclones")
        mkpath(joinpath(outdir, sub))
    end

    ids = dom.reefs.reef_siteid
    write_gpkg(joinpath(outdir, "spatial", "$name.gpkg"), dom.reefs)
    write_cover_nc(joinpath(outdir, "spatial", "coral_cover.nc"), ids, dom.cover)
    write_connectivity(
        joinpath(outdir, "connectivity", connectivity_filename(dom.config)),
        ids, dom.connectivity, dom.config
    )
    for label in climate_labels(dom)
        write_env_nc(
            joinpath(outdir, "DHWs", "dhwRCP$label.nc"), "dhw", dom.dhw[label], dom.reefs;
            long_name="Degree heating weeks (synthetic, climate label $label)",
            units="degC-weeks"
        )
        write_env_nc(
            joinpath(outdir, "waves", "wave_RCP$label.nc"), "Ub", dom.waves, dom.reefs;
            long_name="Wave orbital velocity (synthetic, constant per reef)", units="m s-1"
        )
    end
    write_cyclone_nc(
        joinpath(outdir, "cyclones", "cyclone_mortality.nc"), ids, dom.cyclone_mortality
    )
    write_datapackage(joinpath(outdir, "datapackage.json"), name, dom)
    save_config(joinpath(outdir, "generator_config.toml"), dom.config)
    write_domain_readme(joinpath(outdir, "README.md"), name, dom)

    return outdir
end

"""
    generate_domain(cfg_or_path; outdir=joinpath("domains", <name>), force=false)::String

Build and write a domain in one step, from a `DomainConfig` or the path to a TOML config.
Returns the output folder, ready for `ADRIA.load_domain(outdir, <climate label>)`.
"""
function generate_domain(cfg::DomainConfig; outdir=joinpath("domains", cfg.name), force=false)
    return write_domain(build_domain(cfg), outdir; force=force)
end
function generate_domain(path::AbstractString; kwargs...)
    return generate_domain(load_config(path); kwargs...)
end

reef_outlines(dom::SyntheticDomain; units::Symbol=:lonlat) =
    reef_outlines(dom.reefs, dom.config.reefs.center; units=units)

"""
    climate_labels(domain_path::AbstractString)::Vector{Int}

Climate labels available in a written domain (from its `DHWs/dhwRCP<label>.nc` files), for
use with `ADRIA.load_domain(domain_path, label)`.
"""
function climate_labels(domain_path::AbstractString)::Vector{Int}
    dhw_dir = joinpath(domain_path, "DHWs")
    isdir(dhw_dir) || throw(ArgumentError("No DHWs folder in $domain_path."))
    matches = (match(r"^dhwRCP(\d+)\.nc$", f) for f in readdir(dhw_dir))
    return sort!([parse(Int, m.captures[1]) for m in matches if m !== nothing])
end
