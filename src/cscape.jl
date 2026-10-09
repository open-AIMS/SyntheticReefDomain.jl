# Export a synthetic domain as an input project for the original R package cscape
# (https://github.com/open-AIMS/c-scape), which runs with `cscape::run_cscape(id, fpath)`.
#
# Julia cannot write R's .rds format, so the export has two steps:
#   1. Julia writes the inputs as CSV (inputs/), the initial cover CSV C~scape reads directly
#      and an R script, build_cscape_project.R.
#   2. The R script turns inputs/ into ScenarioID.xlsx and data/*.rds using cscape's own
#      helpers (create_spatial_data, create_scenario_table) and saveRDS. It runs
#      automatically when Rscript is available; otherwise run it later on a machine with R.
# The script *calls* cscape (GPL-3); no cscape code is copied into this package.
#
# The translation follows CscapeInterface.jl/src/domain_translator.jl (the in-memory
# ADRIA -> C~scape.jl translator) so both routes agree:
#   - site table: reef_siteid, k (%), area (m2), depth_med, ub_med (+ long, lat)
#   - connectivity: transpose of our [source, sink] matrix for the 5 broadcast spawners;
#     a local distance-based matrix for the brooder
#   - DHW: per reef in slot 2 (DHW_site); cyclones as categories 0-5 in slot 3
#   - initial cover: % of site area; corymbose non-Acropora split evenly with the brooder
# Demography files (demog_inputs_<type>.rds, fitted brms models) are provided separately.

"""C~scape's functional types, in the order its files and scenario table use."""
const CSCAPE_TYPES = [
    "acro_table", "acro_corym", "corym_non_acro", "small_massive", "large_massive", "brooder"
]
"""ADRIA group (index into `FUNCTIONAL_GROUPS`) used for each C~scape type."""
const CSCAPE_TYPE_TO_GROUP = [1, 2, 3, 4, 5, 3]
const CSCAPE_DISTURBANCES = [
    "DHW", "DHW_site", "cyclone", "COTS_bin1", "COTS_bin2", "COTS_bin3", "COTS_bin4"
]
"""Regions (`Cluster` levels) of the current C~scape demography models."""
const CSCAPE_REGIONS = ["Offshore_Northern", "Offshore_Central", "Offshore_Southern"]
# Covariate ranges the C~scape demography models were fitted on
const CSCAPE_DEPTH_RANGE = (2.9, 12.0)
const CSCAPE_UB_RANGE = (0.047, 0.378)

"""
    CscapeOptions(; region, demography_files, connectivity_scale, brood_site_value, brooder_neighbour_km)

Settings for exporting a domain to C~scape (R). Kept separate from `DomainConfig`: the same
domain can be exported with different options.

- `region` : C~scape `Region`, which must be a `Cluster` level of the demography models
  (currently $(join(CSCAPE_REGIONS, ", ")))
- `demography_files` : the six demography file names (in `CSCAPE_TYPES` order) that you will
  place in the project's `data/` folder; written to the scenario table's `Growth_Surv_file`
- `connectivity_scale` : multiplier for the spawner connectivity. C~scape uses absolute
  transfer probabilities (ADRIA only uses proportions), so you may want to scale ours down
- `brood_site_value` : total dispersal probability of the brooder per source reef (90% stays
  home, 10% goes to reefs within `brooder_neighbour_km`)
- `brooder_neighbour_km` : neighbourhood radius for brooder dispersal
"""
Base.@kwdef struct CscapeOptions
    region::String = "Offshore_Central"
    demography_files::Vector{String} = ["demog_inputs_$(t).rds" for t in CSCAPE_TYPES]
    connectivity_scale::Float64 = 1.0
    brood_site_value::Float64 = 0.03
    brooder_neighbour_km::Float64 = 1.0
end

function validate(opts::CscapeOptions)::CscapeOptions
    errs = String[]
    isempty(strip(opts.region)) && push!(errs, "region must not be empty.")
    if length(opts.demography_files) != length(CSCAPE_TYPES)
        push!(errs, "demography_files needs $(length(CSCAPE_TYPES)) names " *
                    "($(join(CSCAPE_TYPES, ", "))).")
    elseif any(f -> isempty(strip(f)) || occursin(r"[/\\]", f), opts.demography_files)
        push!(errs, "demography_files must be plain file names (no folders): they are " *
                    "looked up in the project's data/ folder.")
    end
    opts.connectivity_scale >= 0 || push!(errs, "connectivity_scale must be >= 0.")
    0 <= opts.brood_site_value <= 1 || push!(errs, "brood_site_value must lie in [0, 1].")
    opts.brooder_neighbour_km >= 0 || push!(errs, "brooder_neighbour_km must be >= 0.")
    isempty(errs) || throw(ArgumentError(
        "Invalid C~scape export options:\n  - " * join(errs, "\n  - ")
    ))
    return opts
end

cscape_options_to_dict(o::CscapeOptions) = _section_dict(o)

"""
    load_cscape_options(path_or_io_or_dict)::CscapeOptions

Read C~scape export options from TOML (e.g. a project's `cscape_export.toml`).
"""
load_cscape_options(path::AbstractString) = load_cscape_options(TOML.parsefile(path))
load_cscape_options(io::IO) = load_cscape_options(TOML.parse(io))
function load_cscape_options(d::AbstractDict)::CscapeOptions
    section = get(d, "cscape", d)
    return validate(_from_dict(CscapeOptions, section, "[cscape]"))
end

"""Write C~scape export options as TOML (`[cscape]` section)."""
function save_cscape_options(path::AbstractString, opts::CscapeOptions)
    open(path, "w") do io
        println(io, "# C~scape export options (SyntheticReefDomain.jl)")
        TOML.print(io, Dict("cscape" => cscape_options_to_dict(opts)); sorted=true)
    end
    return path
end

# --------------------------------------------------------------------------------------
# Translation
# --------------------------------------------------------------------------------------

"""
    cyclone_category(intensity)::Int

C~scape cyclone category for an event intensity in (0, 1]: `ceil(5 x intensity)`, clamped to
1-5 (0.2 -> 1, ..., 1.0 -> 5); 0 when there is no cyclone.
"""
function cyclone_category(intensity::Real)::Int
    intensity > 0 || return 0
    return clamp(ceil(Int, round(5 * intensity; digits=9)), 1, 5)
end

"""
    cyclone_scenario_for_replicate(r, n_cyclone_files)

Cyclone scenario used for C~scape replicate `r` (= DHW scenario `r`). Replicate 1 is the
disturbance-free control; when there are fewer cyclone than DHW scenarios, the stochastic
cyclone scenarios are reused in order.
"""
function cyclone_scenario_for_replicate(r::Int, n_cyclone_files::Int)::Int
    r <= n_cyclone_files && return r
    return ZERO_SCENARIO + 1 + mod(r - ZERO_SCENARIO - 1, n_cyclone_files - ZERO_SCENARIO)
end

"""Brooder dispersal [destination, source]: 90% retained, 10% to neighbours within `radius_km`."""
function brooder_connectivity(reefs::DataFrame, value::Float64, radius_km::Float64)::Matrix{Float64}
    n = nrow(reefs)
    dist = distance_matrix_km(reefs)
    brood = zeros(n, n)
    for src in 1:n
        brood[src, src] = 0.9 * value
        neighbours = [d for d in 1:n if d != src && dist[d, src] <= radius_km]
        isempty(neighbours) || (brood[neighbours, src] .= 0.1 * value / length(neighbours))
    end
    return brood
end

"""
    cscape_inputs(dom::SyntheticDomain, opts::CscapeOptions=CscapeOptions())

The C~scape inputs of a domain, in memory (for previews and tests):

- `sites` : site table (`reef_siteid`, `long`, `lat`, `k` %, `area` m2, `depth_med`,
  `ub_med`, `habitat`, `cluster_id`)
- `connectivity` : `[destination, source, type]` array (`n x n x 6`)
- `disturbances` : climate label => `[year, site, replicate, 7]` array (slot 2 DHW_site,
  slot 3 cyclone category; replicate 1 is the disturbance-free control)
- `init_cover` : initial cover table, % of site area per C~scape type
- `cyclone_scenarios` : cyclone scenario used for each replicate
- `warnings` : reefs outside the demography models' covariate ranges
"""
function cscape_inputs(dom::SyntheticDomain, opts::CscapeOptions=CscapeOptions())
    validate(opts)
    reefs = dom.reefs
    n = nrow(reefs)
    ids = collect(String, reefs.reef_siteid)

    sites = DataFrame(
        reef_siteid=ids, long=reefs.lon, lat=reefs.lat, k=reefs.k, area=reefs.area,
        depth_med=reefs.depth_med, ub_med=reefs.ub_med, habitat=reefs.habitat,
        cluster_id=reefs.cluster_id
    )

    spawners = permutedims(dom.connectivity) .* opts.connectivity_scale  # [dest, src]
    connectivity = Array{Float64}(undef, n, n, length(CSCAPE_TYPES))
    for t in 1:(length(CSCAPE_TYPES) - 1)
        connectivity[:, :, t] .= spawners
    end
    connectivity[:, :, end] .= brooder_connectivity(
        reefs, opts.brood_site_value, opts.brooder_neighbour_km
    )

    n_rep = dom.config.n_dhw_scenarios + 1
    n_cyc = size(dom.cyclone_intensity, 2)
    cyclone_scenarios = [cyclone_scenario_for_replicate(r, n_cyc) for r in 1:n_rep]
    categories = cyclone_category.(dom.cyclone_intensity[:, cyclone_scenarios])  # T x reps
    T = length(years(dom))
    disturbances = Dict{Int,Array{Float64,4}}()
    for label in climate_labels(dom)
        A = zeros(T, n, n_rep, length(CSCAPE_DISTURBANCES))
        for r in 1:n_rep, s in 1:n
            A[:, s, r, 1] .= dom.dhw_regional[label][:, r]
            A[:, s, r, 2] .= dom.dhw[label][:, s, r]
            A[:, s, r, 3] .= categories[:, r]
        end
        disturbances[label] = A
    end

    k_pct = reefs.k
    init_cover = DataFrame(reef_siteid=ids)
    for (t, type) in enumerate(CSCAPE_TYPES)
        group = CSCAPE_TYPE_TO_GROUP[t]
        n_sharing = count(==(group), CSCAPE_TYPE_TO_GROUP)
        init_cover[!, type] = dom.cover[group, :] .* k_pct ./ n_sharing
    end

    warnings = String[]
    for (label, values, (lo, hi), unit) in (("depth_med", reefs.depth_med, CSCAPE_DEPTH_RANGE, "m"),
        ("ub_med", reefs.ub_med, CSCAPE_UB_RANGE, "m/s"))
        outside = findall(v -> !(lo <= v <= hi), values)
        isempty(outside) || push!(warnings,
            "$label outside the demography models' range ($lo-$hi $unit) on " *
            join(["$(ids[i]) ($(round(values[i]; digits=3)))" for i in outside], ", ") *
            ". C~scape will extrapolate for these reefs."
        )
    end

    return (; sites, connectivity, disturbances, init_cover, cyclone_scenarios, warnings)
end

# --------------------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------------------

_csv_field(x::AbstractString) = "\"" * replace(x, "\"" => "\"\"") * "\""
_csv_field(x::Integer) = string(x)
_csv_field(x::Real) = string(Float64(x))
_csv_field(::Nothing) = ""

function _write_csv(path::AbstractString, header::Vector{String}, rows)
    open(path, "w") do io
        println(io, join(_csv_field.(header), ","))
        for row in rows
            println(io, join(_csv_field.(row), ","))
        end
    end
    return path
end

_cscape_files(name, label) = (
    spatial="reef_sites_$(name).rds",
    connectivity="connectivity_$(name).rds",
    temporal="temporal_$(name)_$(label).rds"
)

# Columns of the scenario table, following cscape::create_example_project(), plus `RCP`
# for the older C_scape scripts. Values are written as text and typed by the R script.
function _scenario_rows(dom::SyntheticDomain, name::AbstractString, folder::AbstractString, opts::CscapeOptions)
    cfg = dom.config
    n_rep = cfg.n_dhw_scenarios + 1
    header = [
        "ID", "folder", "num_repeats", "Purpose", "Intervention", "date", "In date", "Region",
        "Year_Start", "Year_End", "init_cover", "fts", "temp_growth", "HeatToleranceGroups",
        "HeatToleranceInit", "Heritability", "Plasticity", "TradeOff", "counterfactual",
        "Spatial_file", "Connectivity_file", "Disturbance_file", "Geometry_file",
        "Growth_Surv_file", "Fogging_reducer", "Fogging_sites", "Fogging_start",
        "Deployment area", "TotalCorals", "species", "species_proportions", "Enhancement",
        "InterventionYears_start", "duration", "frequency", "Reef_siteids", "RCP"
    ]
    rows = map(climate_labels(cfg)) do label
        f = _cscape_files(name, label)
        values = Dict{String,Any}(
            "ID" => label,
            "folder" => folder,
            "num_repeats" => "1:$n_rep",
            "Purpose" => "Synthetic domain $name (SyntheticReefDomain.jl), climate label $label",
            "Intervention" => "No",
            "date" => Dates.format(today(), "dd/mm/yyyy"),
            "In date" => "Yes",
            "Region" => opts.region,
            "Year_Start" => cfg.timeframe[1],
            "Year_End" => cfg.timeframe[2],
            "init_cover" => "initialisationData.csv",
            "fts" => join(CSCAPE_TYPES, "/"),
            "temp_growth" => "FALSE",
            "HeatToleranceGroups" => "-1_1_3",
            "HeatToleranceInit" => "0_0",
            "Heritability" => "1_0.01",
            "Plasticity" => 0,
            "Spatial_file" => f.spatial,
            "Connectivity_file" => f.connectivity,
            "Disturbance_file" => f.temporal,
            "Growth_Surv_file" => join(opts.demography_files, "/"),
            "RCP" => label
        )
        [get(values, h, nothing) for h in header]
    end
    return header, rows
end

const CSCAPE_BUILD_SCRIPT = raw"""
# Builds the C~scape project in this folder from the CSV files in inputs/.
# Generated by SyntheticReefDomain.jl. Requires the R package cscape
# (https://github.com/open-AIMS/c-scape) and its dependencies (sf, writexl, ...).
#
#   Rscript build_cscape_project.R      (works from any folder)
#
# Writes data/reef_sites_*.rds, data/connectivity_*.rds, data/temporal_*.rds and
# ScenarioID.xlsx. data/initialisationData.csv is already in place. Demography files
# (Growth_Surv_file in the scenario table) are not created here: add them to data/.

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grep("^--file=", args)])
root <- if (length(script) == 1) dirname(normalizePath(script)) else getwd()
inp <- function(f) file.path(root, "inputs", f)
dat <- function(f) file.path(root, "data", f)

for (pkg in c("cscape", "sf", "writexl")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("R package '", pkg, "' is required (install cscape from ",
         "https://github.com/open-AIMS/c-scape).", call. = FALSE)
  }
}
for (d in c("data", "logs", "model_outputs", "plots")) {
  dir.create(file.path(root, d), showWarnings = FALSE, recursive = TRUE)
}

scen <- read.csv(inp("scenarios.csv"), colClasses = "character", check.names = FALSE)
scen[] <- lapply(scen, function(x) { x[x == ""] <- NA; x })
sites <- read.csv(inp("sites.csv"), stringsAsFactors = FALSE)
sites$reef_siteid <- as.character(sites$reef_siteid)
ids <- sites$reef_siteid
fts <- strsplit(scen$fts[1], "/", fixed = TRUE)[[1]]
dist_types <- c("DHW", "DHW_site", "cyclone",
                "COTS_bin1", "COTS_bin2", "COTS_bin3", "COTS_bin4")

# 1. Site table (sf points, as cscape writes it)
cscape::create_spatial_data(
  data = sites, coords = c("long", "lat"), output_formats = "rds",
  output_file = sub("\\.rds$", "", scen$Spatial_file[1]),
  output_dir = file.path(root, "data")
)

# 2. Connectivity: [destination, source, coral type]
cl <- read.csv(inp("connectivity.csv"), stringsAsFactors = FALSE)
conn <- array(0, dim = c(length(ids), length(ids), length(fts)),
              dimnames = list(ids, ids, fts))
conn[cbind(match(cl$destination, ids), match(cl$source, ids), match(cl$type, fts))] <- cl$value
saveRDS(conn, dat(scen$Connectivity_file[1]))

# 3. Temporal files: one per climate label (scenario row)
for (i in seq_len(nrow(scen))) {
  dl <- read.csv(inp(paste0("disturbances_", scen$RCP[i], ".csv")), stringsAsFactors = FALSE)
  years <- as.character(sort(unique(dl$year)))
  reps <- sort(unique(dl$replicate))
  dn <- function(vars, var_name) {
    stats::setNames(list(years, ids, paste0("rep_", reps), vars),
                    c("year", "site", "replicate", var_name))
  }
  dims <- c(length(years), length(ids), length(reps))
  disturbances <- array(0, dim = c(dims, length(dist_types)),
                        dimnames = dn(dist_types, "disturbance"))
  idx <- cbind(match(as.character(dl$year), years), match(dl$site, ids),
               match(dl$replicate, reps))
  disturbances[cbind(idx, 1)] <- dl$DHW
  disturbances[cbind(idx, 2)] <- dl$DHW_site
  disturbances[cbind(idx, 3)] <- dl$cyclone
  temporal <- list(
    disturbances = disturbances,
    eggs = array(0, dim = c(dims, length(fts)), dimnames = dn(fts, "coral_type")),
    temp_growth = array(1, dim = c(dims, length(fts)), dimnames = dn(fts, "coral_type"))
  )
  saveRDS(temporal, dat(scen$Disturbance_file[i]))
}

# 4. Scenario table, with the cell types C~scape expects
for (col in c("ID", "Year_Start", "Year_End", "Plasticity")) {
  scen[[col]] <- as.numeric(scen[[col]])
}
scen$temp_growth <- as.logical(scen$temp_growth)
cscape::create_scenario_table(scen, file.path(root, "ScenarioID.xlsx"))

missing <- setdiff(strsplit(scen$Growth_Surv_file[1], "/", fixed = TRUE)[[1]],
                   list.files(file.path(root, "data")))
message("C~scape project built in ", root)
if (length(missing) > 0) {
  message("Add these demography files to ", file.path(root, "data"), " before running:\n  ",
          paste(missing, collapse = "\n  "))
}
"""

const _RSCRIPT_STATUS = Ref{Any}(nothing)

"""
    rscript_status(; refresh=false)

`(rscript, cscape)`: the path to `Rscript` (or `nothing`) and whether the R package cscape is
installed. Cached after the first call.
"""
function rscript_status(; refresh::Bool=false)
    (_RSCRIPT_STATUS[] === nothing || refresh) || return _RSCRIPT_STATUS[]
    rscript = Sys.which("Rscript")
    if rscript === nothing && haskey(ENV, "R_HOME")
        candidate = joinpath(ENV["R_HOME"], "bin", Sys.iswindows() ? "Rscript.exe" : "Rscript")
        isfile(candidate) && (rscript = candidate)
    end
    has_cscape = false
    if rscript !== nothing
        has_cscape = try
            strip(read(`$rscript -e "cat(requireNamespace('cscape', quietly = TRUE))"`, String)) == "TRUE"
        catch
            false
        end
    end
    _RSCRIPT_STATUS[] = (rscript=rscript, cscape=has_cscape)
    return _RSCRIPT_STATUS[]
end

"""
    write_cscape_project(dom, dir=joinpath("cscape_projects", dom.config.name); options=CscapeOptions(), run_r=true, force=false, figures=true)

Write `dom` as an input project for the R package cscape (see the file header for the
layout and translation). With `run_r=true` and `Rscript` + cscape available, the generated
`build_cscape_project.R` is run immediately; otherwise only the CSV inputs and the script
are written. Demography files are never written: add them to `<dir>/data/`.

The project's `README.md` says how to run it and reports every domain setting, a per-reef
table and, when CairoMakie is loaded and `figures=true`, the domain and C~scape plots (saved
in `figures/`).

An existing export is only updated with `force=true`; files you added (e.g. demography) are
kept. Returns `(path, built, message, warnings)`.
"""
function write_cscape_project(
    dom::SyntheticDomain, dir::AbstractString=joinpath("cscape_projects", dom.config.name);
    options::CscapeOptions=CscapeOptions(), run_r::Bool=true, force::Bool=false,
    figures::Bool=true
)
    validate(options)
    dir = normpath(abspath(dir))
    name = dom.config.name

    if isdir(dir) && !isempty(readdir(dir))
        is_export = isfile(joinpath(dir, "build_cscape_project.R"))
        is_export || throw(ArgumentError(
            "Refusing to write into $dir: it is not empty and is not a C~scape export of this " *
            "package (no build_cscape_project.R)."
        ))
        force || throw(ArgumentError(
            "$dir already holds a C~scape export. Pass force=true to update it (files you " *
            "added, such as demography files, are kept)."
        ))
        # Remove only files this function writes, so stale climate labels don't linger
        for (sub, pattern) in (("inputs", r"\.csv$"), ("data", Regex("^temporal_\\Q$(name)\\E_\\d+\\.rds\$")))
            folder = joinpath(dir, sub)
            isdir(folder) || continue
            foreach(f -> occursin(pattern, f) && rm(joinpath(folder, f)), readdir(folder))
        end
        rm(joinpath(dir, "figures"); recursive=true, force=true)
    end

    inputs = cscape_inputs(dom, options)
    ids = inputs.sites.reef_siteid
    for sub in ("inputs", "data", "logs", "model_outputs", "plots")
        mkpath(joinpath(dir, sub))
    end

    s = inputs.sites
    _write_csv(joinpath(dir, "inputs", "sites.csv"), names(s),
        (collect(row) for row in eachrow(s)))

    C = inputs.connectivity
    _write_csv(joinpath(dir, "inputs", "connectivity.csv"),
        ["destination", "source", "type", "value"],
        ((ids[d], ids[src], CSCAPE_TYPES[t], C[d, src, t])
         for t in axes(C, 3), src in axes(C, 2), d in axes(C, 1) if C[d, src, t] != 0))

    yrs = collect(years(dom))
    for (label, A) in inputs.disturbances
        _write_csv(joinpath(dir, "inputs", "disturbances_$(label).csv"),
            ["year", "site", "replicate", "DHW", "DHW_site", "cyclone"],
            ((yrs[t], ids[s], r, A[t, s, r, 1], A[t, s, r, 2], round(Int, A[t, s, r, 3]))
             for r in axes(A, 3), s in axes(A, 2), t in axes(A, 1)))
    end

    header, rows = _scenario_rows(dom, name, basename(dir), options)
    _write_csv(joinpath(dir, "inputs", "scenarios.csv"), header, rows)

    ic = inputs.init_cover
    _write_csv(joinpath(dir, "data", "initialisationData.csv"), names(ic),
        (collect(row) for row in eachrow(ic)))

    write(joinpath(dir, "build_cscape_project.R"), CSCAPE_BUILD_SCRIPT)
    save_cscape_options(joinpath(dir, "cscape_export.toml"), options)
    save_config(joinpath(dir, "generator_config.toml"), dom.config)

    built, message = false, ""
    status = run_r ? rscript_status() : (rscript=nothing, cscape=false)
    if !run_r
        message = "Inputs and build script written (R step skipped)."
    elseif status.rscript === nothing
        message = "Rscript not found: inputs and build script written. Run " *
                  "`Rscript build_cscape_project.R` on a machine with R and cscape."
    elseif !status.cscape
        message = "R found, but the cscape package is not installed: inputs and build " *
                  "script written. Install cscape, then run `Rscript build_cscape_project.R`."
    else
        log = IOBuffer()
        script = joinpath(dir, "build_cscape_project.R")
        try
            run(pipeline(Cmd(`$(status.rscript) $script`; dir=dir); stdout=log, stderr=log))
            built = true
            message = "C~scape project built with R."
        catch err
            message = "The R build step failed (see build_log.txt): " * sprint(showerror, err)
        end
        write(joinpath(dir, "build_log.txt"), take!(log))
    end

    figure_files = if !figures
        nothing
    elseif plots_available()
        # invokelatest: the extension may have been loaded after this code was compiled
        Base.invokelatest(save_cscape_figures, dom, options, joinpath(dir, "figures"))
    else
        @info "CairoMakie is not loaded, so the C~scape project's README has no plots. Run " *
              "`using CairoMakie` before writing to include them."
        nothing
    end
    _write_cscape_readme(joinpath(dir, "README.md"), dom, options, inputs, built, message, figure_files)
    return (path=dir, built=built, message=message, warnings=inputs.warnings)
end

function _write_cscape_readme(path, dom::SyntheticDomain, opts::CscapeOptions, inputs, built::Bool,
                              message::AbstractString, figures)
    cfg = dom.config
    labels = climate_labels(cfg)
    n_rep = cfg.n_dhw_scenarios + 1
    fpath = replace(dirname(path), '\\' => '/')
    summary = sprint(show, MIME("text/plain"), dom)
    open(path, "w") do io
        println(io, "# C~scape project: $(cfg.name)\n")
        println(io, "Input project for the R package [cscape](https://github.com/open-AIMS/c-scape), " *
                    "generated by $GENERATOR on $(Dates.format(now(), "yyyy-mm-dd HH:MM")) from " *
                    "the synthetic domain `$(cfg.name)` (settings in `generator_config.toml`, " *
                    "export options in `cscape_export.toml`). Not observational data: see the " *
                    "package README's *Provenance* section.\n")
        println(io, "The matching ADRIA domain comes from the same settings: run " *
                    "`generate_domain(\"generator_config.toml\")` with SyntheticReefDomain.jl, or " *
                    "open that file in the `reef_builder.jl` notebook and write the domain.\n")
        println(io, "```text\n$summary\n```\n")
        println(io, "**Status:** $message\n")
        built || println(io, "Build the R files first: `Rscript build_cscape_project.R`.\n")
        println(io, "## Before running\n")
        println(io, "Add the six demography files to `data/` (they are not part of this export):\n")
        for (t, f) in zip(CSCAPE_TYPES, opts.demography_files)
            println(io, "- `$f` ($t)")
        end
        println(io, "\nRegion (must be a `Cluster` level of those models): `$(opts.region)`.\n")
        println(io, "## Running\n")
        println(io, "```r\nlibrary(cscape)")
        for label in labels
            println(io, "run_cscape($label, \"$fpath\", rep_index = 2)   # one replicate, climate label $label")
        end
        println(io, "run_cscape_batch($(first(labels)), \"$fpath\")       # all $n_rep replicates\n```\n")
        println(io, "Scenario IDs are the climate labels. Replicate r uses DHW scenario r of " *
                    "the domain; **replicate 1 is the disturbance-free control**.\n")
        println(io, "## How the domain was translated\n")
        println(io, "- Sites: `reef_siteid`, `k` (%), `area` (m²), `depth_med`, `ub_med`, plus `long`/`lat`, `habitat`, `cluster_id`.")
        println(io, "- Connectivity: `[destination, source, type]`; the 5 spawning types use the domain's matrix transposed " *
                    "(x $(opts.connectivity_scale)); the brooder keeps $(round(90 * opts.brood_site_value; digits=2))% at home and " *
                    "spreads $(round(10 * opts.brood_site_value; digits=2))% over reefs within $(opts.brooder_neighbour_km) km. " *
                    "C~scape uses these values as absolute transfer probabilities.")
        println(io, "- Heat stress: per-reef DHW in slot 2 (`DHW_site`); the regional signal in slot 1 (unused by C~scape).")
        println(io, "- Cyclones: categories 1-5 from the event intensity (`ceil(5 x intensity)`), the same on every reef. " *
                    "C~scape applies its own wind-speed mortality, so the domain's per-group cyclone vulnerability is not used.")
        println(io, "- Initial cover: % of site area; corymbose non-Acropora is split evenly between `corym_non_acro` and `brooder`.")
        println(io, "- No COTS, no external eggs (closed system), `temp_growth` = 1.\n")
        if !isempty(inputs.warnings)
            println(io, "## Warnings\n")
            foreach(w -> println(io, "- $w"), inputs.warnings)
            println(io)
        end

        # The domain report: every setting, the per-reef table and the plots
        _settings(io, cfg.name, dom; adria=false)
        println(io, "### C~scape export options\n")
        _table(io, ["Setting", "Value"], [
            ["Region", "`$(opts.region)`"],
            ["Demography files", join(["`$f`" for f in opts.demography_files], ", ")],
            ["Spawner connectivity scale", _num(opts.connectivity_scale)],
            ["Brooder dispersal per source reef", "$(_num(opts.brood_site_value)) " *
                "(90% stays home, 10% to reefs within $(_num(opts.brooder_neighbour_km)) km)"]
        ])
        _plots(io, figures)
        _cscape_files_table(io, cfg, figures)
    end
    return path
end

function _cscape_files_table(io::IO, cfg::DomainConfig, figures)
    println(io, "## Files\n")
    labels = climate_labels(cfg)
    f = _cscape_files(cfg.name, first(labels))
    rows = Any[
        ["`ScenarioID.xlsx`", "Scenario table, one row per climate label (built by the R script)"],
        ["`data/$(f.spatial)`", "Site table (built by the R script)"],
        ["`data/$(f.connectivity)`", "Connectivity `[destination, source, type]` (built by the R script)"]
    ]
    for label in labels
        push!(rows, ["`data/$(_cscape_files(cfg.name, label).temporal)`",
            "Disturbances `[year, site, replicate, 7]`, climate label $label (built by the R script)"])
    end
    append!(rows, Any[
        ["`data/initialisationData.csv`", "Initial cover per site and type (% of site area)"],
        ["`data/demog_inputs_<type>.rds`", "Demography models: **not included**, add them yourself"],
        ["`inputs/*.csv`", "The inputs as CSV, read by `build_cscape_project.R`"],
        ["`build_cscape_project.R`", "Builds `ScenarioID.xlsx` and the `.rds` files with cscape's helpers"],
        ["`generator_config.toml`", "The exact settings of the synthetic domain"],
        ["`cscape_export.toml`", "The C~scape export options"],
        ["`logs/`, `model_outputs/`, `plots/`", "Empty folders that C~scape writes into"]
    ])
    figures === nothing || push!(rows, ["`figures/`", "The plots shown above (PNG)"])
    _table(io, ["Path", "Content"], rows)
    return nothing
end

# --------------------------------------------------------------------------------------
# Zip downloads
# --------------------------------------------------------------------------------------

"""Zip the contents of `dir` under the folder name `prefix`; returns the zip file's bytes."""
function zip_folder(dir::AbstractString, prefix::AbstractString)::Vector{UInt8}
    io = IOBuffer()
    ZipArchives.ZipWriter(io) do w
        for (root, _, files) in walkdir(dir), f in sort(files)
            full = joinpath(root, f)
            rel = replace(relpath(full, dir), '\\' => '/')
            ZipArchives.zip_newfile(w, "$prefix/$rel"; compress=true)
            write(w, read(full))
        end
    end
    return take!(io)
end

"""
    domain_zip(dom; figures=true)::Vector{UInt8}

The domain as a zip file (folder `<name>/`), e.g. for a download button.
"""
function domain_zip(dom::SyntheticDomain; figures::Bool=true)::Vector{UInt8}
    tmp = mktempdir()
    path = write_domain(dom, joinpath(tmp, dom.config.name); figures=figures)
    return zip_folder(path, dom.config.name)
end

"""
    cscape_zip(dom, options=CscapeOptions(); run_r=true, figures=true)::Vector{UInt8}

The C~scape export of `dom` as a zip file (folder `<name>_cscape/`).
"""
function cscape_zip(dom::SyntheticDomain, options::CscapeOptions=CscapeOptions();
                    run_r::Bool=true, figures::Bool=true)
    tmp = mktempdir()
    folder = "$(dom.config.name)_cscape"
    result = write_cscape_project(dom, joinpath(tmp, folder); options=options, run_r=run_r,
                                  figures=figures)
    return zip_folder(result.path, folder)
end
