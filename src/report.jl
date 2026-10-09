# The domain's README.md: a report of every setting (as in the notebook's controls), a
# per-reef table, the plots (when CairoMakie is loaded) and the files in the package.

const FIGURE_TITLES = Dict(
    "map.png" => ("Reefs and connectivity",
        "Reefs coloured by initial cover and sized by area; outline colour = cluster. " *
        "Arrows show connectivity from source to sink reef, thicker for stronger links " *
        "(links below 5% of the strongest are not drawn)."),
    "connectivity.png" => ("Connectivity matrix",
        "Transfer probabilities from source reefs (rows) to receiving reefs (columns)."),
    "cover.png" => ("Initial cover", "Initial cover per reef by functional group."),
    "cyclones.png" => ("Cyclones",
        "Left: cyclone years and intensities in scenario 2 (the first stochastic scenario). " *
        "Right: mortality per functional group for an average and the strongest event; " *
        "every reef loses the same fraction."),
    "cscape_connectivity.png" => ("C~scape connectivity",
        "Transfer probabilities as written for C~scape (destination rows, source columns). " *
        "The five broadcast-spawning types share the domain's matrix (times the connectivity " *
        "scale); the brooder keeps 90% of its dispersal at home."),
    "cscape_cyclones.png" => ("C~scape cyclone categories",
        "Cyclone category (0-5) per year and replicate, the same on every reef. Replicate 1 " *
        "is the disturbance-free control."),
    "cscape_cover.png" => ("C~scape initial cover",
        "Initial cover per reef for the six C~scape types, in % of site area " *
        "(`data/initialisationData.csv`).")
)
const DHW_CAPTION =
    "Left: DHW of every reef (colour = exposure) and the regional signal in scenario 2 " *
    "(the first stochastic scenario). Right: mean DHW over reefs across all stochastic " *
    "scenarios (median and 10-90% band) for every climate level."

_num(x::Integer) = string(x)
function _num(x::Real)
    r = round(x; sigdigits=4)
    return r == round(r) && abs(r) < 1e9 ? string(Int(r)) : string(r)
end
_fixed(x::Real, digits::Int) = Printf.format(Printf.Format("%.$(digits)f"), x)
_nums(v) = join(_num.(v), ", ")
_years_text(v) = isempty(v) ? "none" : join(v, ", ")
_yes(b::Bool) = b ? "yes" : "no"

"""Readable description of a per-reef attribute specification."""
function describe_spec(spec)::String
    s = attr_spec(spec)
    s isa Float64 && return "$(_num(s)) on every reef"
    s isa Vector{Float64} && return "per reef: $(_nums(s))"
    haskey(s, :min) && return "$(_num(s.min)) to $(_num(s.max)), random per reef"
    return "gradient from $(_num(s.from)) (first reef) to $(_num(s.to)) (last reef)"
end
_describe_categorical(spec) =
    spec isa AbstractVector ? "per reef: " * join(string.(spec), ", ") : string(spec)

function _table(io::IO, header, rows)
    println(io, "| ", join(header, " | "), " |")
    println(io, "|", join(fill(" --- ", length(header)), "|"), "|")
    for row in rows
        println(io, "| ", join(string.(row), " | "), " |")
    end
    println(io)
    return nothing
end

# `adria=false` leaves out the ADRIA file names (for the C~scape project's report)
function _settings(io::IO, name::AbstractString, dom::SyntheticDomain; adria::Bool=true)
    cfg = dom.config
    r, c, cv, d, cy = cfg.reefs, cfg.connectivity, cfg.cover, cfg.dhw, cfg.cyclones
    labels = climate_labels(cfg)
    T = length(years(cfg))
    n = r.n
    stochastic(A) = selectdim(A, ndims(A), (ZERO_SCENARIO + 1):size(A, ndims(A)))

    println(io, "## Settings\n")

    println(io, "### Domain\n")
    _table(io, ["Setting", "Value"], [
        ["Name", name == cfg.name ? name : "$name (configured as $(cfg.name))"],
        ["Years", "$(cfg.timeframe[1])-$(cfg.timeframe[2]) ($T years)"],
        ["Seed", cfg.seed],
        ["Climate labels", adria ? join(["$l (`DHWs/dhwRCP$l.nc`)" for l in labels], ", ") :
                           join(labels, ", ")],
        ["DHW scenarios", "1 control (no heat stress) + $(cfg.n_dhw_scenarios) stochastic"],
        ["Cyclone scenarios", "1 control (no cyclones) + $(cfg.n_cyclone_scenarios) stochastic"]
    ])

    println(io, "### Reefs\n")
    lon, lat = r.center
    rows = Any[
        ["Number of reefs", n],
        ["Layout", "$(r.layout), main axis $(_num(r.bearing_deg))° from north"],
        ["Spacing", "$(_num(r.spacing_km)) km"],
        ["Clusters", r.n_clusters]
    ]
    r.layout == "clustered" &&
        push!(rows, ["Spacing between cluster centres", "$(_num(r.cluster_spacing_km)) km"])
    append!(rows, Any[
        ["Centre", "$(round(abs(lon); digits=3))° $(lon >= 0 ? "E" : "W"), " *
                   "$(round(abs(lat); digits=3))° $(lat >= 0 ? "N" : "S")"],
        ["Area (m²)", describe_spec(r.area_m2)],
        ["Habitable fraction k (%)", describe_spec(r.k_percent)],
        ["Median depth (m)", describe_spec(r.depth_med)],
        ["Median wave orbital velocity `ub_med` (m/s)", describe_spec(r.ub_med)],
        ["Habitat", _describe_categorical(r.habitat)],
        ["CoralBlox calibration group", _describe_categorical(r.cb_calib_group)],
        ["Nudge identical (ub_med, depth_med) pairs apart", _yes(r.distinct_ub_depth)]
    ])
    adria && push!(rows, ["Wave file", r.waves_from_ub_med ? "each reef's `ub_med`, constant" : "zeros"])
    _table(io, ["Setting", "Value"], rows)

    println(io, "### Connectivity\n")
    C = dom.connectivity
    rows = if c.topology == "truncated_dilution"
        Any[
            ["Topology", "truncated_dilution: distance-based connectivity for C~scape " *
                         "(truncated exponential kernel, weighted by destination substrate area)"],
            ["No settlement within", "$(c.threshold_m) m of a source reef's centre"],
            ["Dispersal distance (kernel down to 5%)", "$(_num(c.dispersal_km)) km"],
            ["Matrix total", "$(_num(c.matrix_total)) set; $(_num(sum(C))) after the 1e-6 cutoff"],
            ["Reefs with local retention",
                "$(count(i -> C[i, i] > 0, 1:n)) of $n (needs a substrate radius above $(c.threshold_m) m)"]
        ]
    else
        Any[
            ["Topology", c.topology],
            ["Self-recruitment", _num(c.self_recruitment)],
            ["Strength (largest off-diagonal value)", _num(c.strength)]
        ]
    end
    if c.topology in ("chain", "ring")
        push!(rows, ["Asymmetry (0 symmetric, 1 downstream only)", _num(c.asymmetry)])
        push!(rows, ["Reach (downstream neighbours linked)", c.reach])
    elseif c.topology == "distance_decay"
        push!(rows, ["Decay distance", "$(_num(c.decay_km)) km"])
    elseif c.topology == "random"
        push!(rows, ["Link density", _num(c.density)])
    end
    push!(rows, ["Off-diagonal links", offdiag_links(C)])
    _table(io, ["Setting", "Value"], rows)

    println(io, "### Initial cover\n")
    basis = cv.basis == "habitable" ? "habitable area" : "whole reef area"
    _table(io, ["Setting", "Value"], [
        ["Total cover (%)", "$(describe_spec(cv.total_percent)), as % of $basis"],
        ["Per-reef variation of the mix (jitter)", _num(cv.jitter)]
    ])
    w = cv.composition ./ sum(cv.composition)
    realised = vec(mean(dom.cover; dims=2)) .* 100
    _table(io, ["Functional group", "Weight", "Share (%)", "Mean initial cover (% of habitable area)"],
        [[FUNCTIONAL_GROUPS[g], _num(cv.composition[g]), _fixed(100w[g], 1), _fixed(realised[g], 2)]
         for g in 1:N_GROUPS])

    println(io, "### Heat stress (DHW)\n")
    _table(io, ["Setting", "Value"], [
        ["Baseline (DHW in the first year)", _num(d.baseline)],
        ["Year-to-year variability (SD)", _num(d.interannual_sd)],
        ["Autocorrelation (AR1)", _num(d.ar1)],
        ["Marine heatwave probability per year", _num(d.event_prob)],
        ["Marine heatwave size (extra DHW)", "$(_num(d.event_magnitude[1])) to $(_num(d.event_magnitude[2]))"],
        ["Forced heatwave years", _years_text(d.event_years)],
        ["Reef exposure (multiplier)", describe_spec(d.reef_exposure)],
        ["Reef offset (DHW)", describe_spec(d.reef_offset)],
        ["Local noise (SD)", _num(d.local_sd)]
    ])
    println(io, "Per-climate-level values (these override the shared settings above). " *
                "Mean DHW is over all reefs and stochastic scenarios.\n")
    level_rows = map(labels) do label
        p = climate_params(d, label)
        overrides = [k for k in sort(collect(keys(d.climate[label]))) if k != "trend_per_decade"]
        other = isempty(overrides) ? "--" :
            join(["$k = $(repr(d.climate[label][k]))" for k in overrides], "; ")
        D = stochastic(dom.dhw[label])
        [label; adria ? ["`DHWs/dhwRCP$label.nc`"] : String[]; _num(p.trend_per_decade); other;
         _fixed(mean(@view D[1:min(10, T), :, :]), 1); _fixed(mean(@view D[max(1, T - 9):T, :, :]), 1)]
    end
    _table(io, ["Label"; adria ? ["File"] : String[]; "Trend (DHW/decade)"; "Other overrides";
                "Mean DHW, first decade"; "Mean DHW, last decade"], level_rows)

    println(io, "### Cyclones (system-wide)\n")
    events = vec(count(>(0), stochastic(dom.cyclone_intensity); dims=1))
    _table(io, ["Setting", "Value"], [
        ["Enabled", _yes(cy.enabled)],
        ["Probability per year", _num(cy.annual_prob)],
        ["Forced cyclone years", _years_text(cy.event_years)],
        ["Event intensity", "$(_num(cy.intensity[1])) to $(_num(cy.intensity[2]))"],
        ["Cyclones per stochastic scenario", "mean $(_fixed(mean(events), 1)), range $(minimum(events))-$(maximum(events))"]
    ])
    lo, hi = cy.intensity
    _table(io, ["Functional group", "Mortality at intensity 1", "Average event", "Strongest event"],
        [[FUNCTIONAL_GROUPS[g], _num(v), _num(clamp(v * (lo + hi) / 2, 0, 1)), _num(clamp(v * hi, 0, 1))]
         for (g, v) in enumerate(cy.group_vulnerability)])

    println(io, "### Reefs (resolved values)\n")
    reefs = dom.reefs
    first_label = first(labels)
    dhw_mean = vec(mean(stochastic(dom.dhw[first_label]); dims=(1, 3)))
    cover_total = vec(sum(dom.cover; dims=1)) .* 100
    _table(io, ["Reef", "Cluster", "Habitat", "Lon", "Lat", "Area (m²)", "k (%)", "Depth (m)",
                "ub_med (m/s)", "DHW exposure", "DHW offset", "Initial cover (%)",
                "Mean DHW (label $first_label)"],
        [[reefs.reef_siteid[i], reefs.cluster_id[i], reefs.habitat[i],
          _fixed(reefs.lon[i], 4), _fixed(reefs.lat[i], 4), _fixed(reefs.area[i], 0),
          _fixed(reefs.k[i], 1), _fixed(reefs.depth_med[i], 1), _fixed(reefs.ub_med[i], 3),
          _fixed(reefs.dhw_exposure[i], 2), _fixed(reefs.dhw_offset[i], 2),
          _fixed(cover_total[i], 1), _fixed(dhw_mean[i], 2)] for i in 1:n])
    return nothing
end

function _plots(io::IO, figures)
    println(io, "## Plots\n")
    if figures === nothing
        println(io, "Plots were not generated (CairoMakie was not loaded, or `figures=false` " *
                    "was passed). Run `using CairoMakie` before writing (the notebooks do this) " *
                    "to include them.\n")
        return nothing
    end
    for file in figures
        m = match(r"^dhw_(\d+)\.png$", file)
        title, caption = m === nothing ? FIGURE_TITLES[file] :
            ("Heat stress, climate label $(m.captures[1])", DHW_CAPTION)
        println(io, "### $title\n")
        println(io, "![$title](figures/$file)\n")
        println(io, "$caption\n")
    end
    return nothing
end

function _files(io::IO, name::AbstractString, dom::SyntheticDomain, figures)
    cfg = dom.config
    println(io, "## Files\n")
    rows = Any[
        ["`datapackage.json`", "ADRIA data package metadata (version $DPKG_VERSION)"],
        ["`spatial/$name.gpkg`", "Reef polygons (EPSG:4326) and attributes"],
        ["`spatial/coral_cover.nc`", "Initial cover, 5 groups x reefs, fraction of habitable area"],
        ["`connectivity/$(connectivity_filename(cfg))`", "Source (row) x sink (column) transfer probabilities"]
    ]
    for label in climate_labels(cfg)
        push!(rows, ["`DHWs/dhwRCP$label.nc`",
            "DHW, years x reefs x $(cfg.n_dhw_scenarios + 1) scenarios (label $label; scenario 1 = no heat stress)"])
        push!(rows, ["`waves/wave_RCP$label.nc`", "Wave Ub, constant per reef"])
    end
    push!(rows, ["`cyclones/cyclone_mortality.nc`",
        "Cyclone mortality, years x reefs x 5 groups x $(cfg.n_cyclone_scenarios + 1) scenarios (scenario 1 = no cyclones)"])
    push!(rows, ["`generator_config.toml`", "The exact settings used to generate this domain"])
    figures === nothing || push!(rows, ["`figures/`", "The plots shown above (PNG)"])
    _table(io, ["Path", "Content"], rows)
    return nothing
end

"""
    write_domain_report(path, name, dom; figures=nothing)

Write the domain README: summary, all settings, a per-reef table, the plots listed in
`figures` (file names inside `figures/`, as returned by `save_figures`) and the file list.
"""
function write_domain_report(
    path::AbstractString, name::AbstractString, dom::SyntheticDomain; figures=nothing
)
    cfg = dom.config
    labels = climate_labels(cfg)
    summary = sprint(show, MIME("text/plain"), dom)
    open(path, "w") do io
        println(io, "# $name\n")
        println(io, "Synthetic ADRIA domain generated by $GENERATOR on " *
                    "$(Dates.format(now(), "yyyy-mm-dd HH:MM")). Not observational data: see " *
                    "the package README's *Provenance* section for where the default values " *
                    "come from.\n")
        println(io, "The settings are stored in `generator_config.toml`. Regenerate this " *
                    "domain with `generate_domain(\"generator_config.toml\")`.\n")
        println(io, "```text\n$summary\n```\n")
        println(io, "Load in ADRIA:\n")
        println(io, "```julia\ndom = ADRIA.load_domain(\"$name\", \"$(first(labels))\")\n```\n")
        println(io, "Scenario 1 of every DHW, wave and cyclone file is a disturbance-free " *
                    "control (zero DHW and zero cyclone mortality in every year); the stochastic " *
                    "scenarios are 2 to n+1.\n")
        _settings(io, name, dom)
        _plots(io, figures)
        _files(io, name, dom, figures)
    end
    return path
end
