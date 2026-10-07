# Configuration for a synthetic domain: one struct per TOML section, conversion to and from
# TOML dictionaries, and validation.
#
# Per-reef attributes (area, k, depth, wave exposure, DHW exposure, initial cover) accept
# four forms so the same config works for 1 or 20 reefs:
#   - a number              every reef gets that value
#   - a list of length n    one value per reef, in reef order
#   - {min = a, max = b}    values drawn uniformly from [a, b] (seeded)
#   - {from = a, to = b}    a linear gradient from the first to the last reef
# Reef order follows the layout (west -> east along a line for bearing 90, serpentine rows
# for a grid, around the ring, ...), so gradients and chain connectivity follow it too.

"""ADRIA functional groups, in the order ADRIA expects them in cover and cyclone files."""
const FUNCTIONAL_GROUPS = [
    "tabular Acropora",
    "corymbose Acropora",
    "corymbose non-Acropora",
    "Small massives",
    "Large massives"
]
const N_GROUPS = length(FUNCTIONAL_GROUPS)

const MAX_REEFS = 20
# Intervention start-year factors in ADRIA go up to year 25. A minimum of 26 years also
# guarantees n_timesteps != n_reefs, which ADRIA's NetCDF loader relies on to find the
# location dimension (it matches dimensions by size).
const MIN_YEARS = 26
const DPKG_VERSION = "0.8.0"  # ADRIA.COMPAT_DPKG

const LAYOUTS = ("line", "grid", "ring", "random", "clustered")
const TOPOLOGIES = ("isolated", "chain", "ring", "all_to_all", "distance_decay", "random")
const COVER_BASES = ("habitable", "reef")
const CLIMATE_KEYS = (
    "baseline", "trend_per_decade", "interannual_sd", "ar1", "event_prob",
    "event_magnitude", "event_years", "local_sd"
)

"""
Reef layout and per-reef attributes (`[reefs]`).

- `n` : number of reefs (1-20); each reef is one ADRIA location
- `layout` : one of $(join(LAYOUTS, ", "))
- `center` : `[lon, lat]` of the reef system centre (EPSG:4326)
- `spacing_km` : distance between neighbouring reefs
- `bearing_deg` : main axis of the layout, clockwise from north (90 = west -> east)
- `n_clusters` : number of `cluster_id` groups (contiguous blocks of reefs)
- `cluster_spacing_km` : distance between cluster centres (`layout = "clustered"` only)
- `area_m2`, `k_percent`, `depth_med`, `ub_med` : per-reef attributes (see file header)
- `habitat` : a habitat name, or one per reef (informational for ADRIA)
- `cb_calib_group` : CoralBlox calibration group (1-12), one value or one per reef
- `distinct_ub_depth` : nudge `ub_med` by 1e-4 m/s where two reefs share an identical
  `(ub_med, depth_med)` pair, because C~scape collapses such reefs onto one demography point
- `waves_from_ub_med` : write each reef's `ub_med` as its (constant) wave Ub time series;
  `false` writes zeros
"""
Base.@kwdef struct ReefConfig
    n::Int = 10
    layout::String = "line"
    center::Vector{Float64} = [146.0, -18.0]
    spacing_km::Float64 = 5.0
    bearing_deg::Float64 = 90.0
    n_clusters::Int = 1
    cluster_spacing_km::Float64 = 20.0
    area_m2::Any = (min=60_000.0, max=93_000.0)
    k_percent::Any = (min=45.0, max=63.0)
    depth_med::Any = (min=2.0, max=15.0)
    ub_med::Any = (min=0.23, max=0.55)
    habitat::Any = "Slope"
    cb_calib_group::Any = 1
    distinct_ub_depth::Bool = true
    waves_from_ub_med::Bool = true
end

"""
Larval connectivity between reefs (`[connectivity]`), as a source (row) x sink (column)
transfer probability matrix.

- `topology` : one of $(join(TOPOLOGIES, ", "))
- `self_recruitment` : diagonal value (fraction of a reef's larvae that settle at home)
- `strength` : largest off-diagonal value
- `asymmetry` : chain/ring only; 0 = symmetric, 1 = forward only (reef i -> i+1)
- `reach` : chain/ring only; number of downstream neighbours linked (weights halve per step)
- `decay_km` : distance_decay only; e-folding distance of exp(-d / decay_km)
- `density` : random only; probability that a directed link exists

Each row is rescaled where `self_recruitment + exports > 1` so that no reef exports more
larvae than it produces. Values below 1e-6 are zeroed (ADRIA applies the same cutoff).
"""
Base.@kwdef struct ConnectivityConfig
    topology::String = "distance_decay"
    self_recruitment::Float64 = 0.2
    strength::Float64 = 0.05
    asymmetry::Float64 = 0.0
    reach::Int = 1
    decay_km::Float64 = 10.0
    density::Float64 = 0.3
end

"""
Initial coral cover (`[cover]`).

- `total_percent` : total cover per reef, a per-reef attribute (see file header)
- `basis` : `"habitable"` = percent of habitable area (ADRIA's native unit),
  `"reef"` = percent of the whole reef area (converted by dividing by k)
- `composition` : relative share of each functional group, in the order of
  `FUNCTIONAL_GROUPS`; normalised to sum to 1
- `jitter` : per-reef variation of the composition (Dirichlet); 0 = same mix on every reef,
  about 0.05 = subtle, 0.5 = strong
"""
Base.@kwdef struct CoverConfig
    total_percent::Any = (min=10.0, max=40.0)
    basis::String = "habitable"
    composition::Vector{Float64} = [0.26, 0.25, 0.07, 0.19, 0.23]
    jitter::Float64 = 0.0
end

"""
Heat stress as degree heating weeks (`[dhw]`).

Each climate level `label => overrides` writes `DHWs/dhwRCP<label>.nc`. The label is a
positive integer: ADRIA stores it as a number in its scenario table, so it is only a key and
does not need to be a real RCP. Overrides may set any of $(join(CLIMATE_KEYS, ", ")).

For scenario `s`, year `t` and reef `r`:

    R[t, s]    = baseline + trend_per_decade * (t - 1) / 10 + a[t, s] + event[t, s]
    DHW[t,r,s] = max(0, reef_exposure[r] * R[t, s] + reef_offset[r] + local_sd * e[t, r, s])

where `a` is an AR(1) anomaly with standard deviation `interannual_sd` and lag-1
autocorrelation `ar1`, `event` adds a marine heatwave of size `event_magnitude = [min, max]`
with probability `event_prob` per year (and always in `event_years`), and `e` is standard
normal noise. Larger `local_sd` lowers the correlation between reefs; `local_sd = 0` with
`reef_exposure = 1` gives identical DHW on every reef. All climate levels share the same
random draws, so scenario `s` differs between levels only through their parameters.
"""
Base.@kwdef struct DHWConfig
    baseline::Float64 = 3.5
    trend_per_decade::Float64 = 1.2
    interannual_sd::Float64 = 2.0
    ar1::Float64 = 0.3
    event_prob::Float64 = 0.15
    event_magnitude::Vector{Float64} = [2.0, 6.0]
    event_years::Vector{Int} = Int[]
    reef_exposure::Any = (min=0.7, max=1.3)
    reef_offset::Any = 0.0
    local_sd::Float64 = 2.0
    climate::Dict{Int,Dict{String,Any}} = Dict(45 => Dict{String,Any}())
end

"""
System-wide cyclones (`[cyclones]`): in a cyclone year every reef suffers the same
mortality, which differs between functional groups.

- `enabled` : `false` writes an all-zero cyclone file
- `annual_prob` : probability of a cyclone in any year (per scenario)
- `event_years` : years with a cyclone in every scenario
- `intensity` : `[min, max]` of the per-event intensity (0-1), drawn uniformly
- `group_vulnerability` : mortality at intensity 1 for each group in `FUNCTIONAL_GROUPS`

Mortality = `clamp(intensity * group_vulnerability[g], 0, 1)`.
"""
Base.@kwdef struct CycloneConfig
    enabled::Bool = true
    annual_prob::Float64 = 0.1
    event_years::Vector{Int} = Int[]
    intensity::Vector{Float64} = [0.2, 1.0]
    group_vulnerability::Vector{Float64} = [0.6, 0.5, 0.4, 0.2, 0.1]
end

"""
    DomainConfig(; name, timeframe, seed, n_dhw_scenarios, n_cyclone_scenarios,
                   reefs, connectivity, cover, dhw, cyclones)

Complete description of a synthetic domain. Build one with keywords, e.g.

```julia
cfg = DomainConfig(;
    name="Chain5",
    reefs=ReefConfig(; n=5),
    connectivity=ConnectivityConfig(; topology="chain", asymmetry=1.0)
)
```

or load one from TOML with [`load_config`](@ref).
"""
Base.@kwdef struct DomainConfig
    name::String = "SyntheticReefs"
    timeframe::Vector{Int} = [2025, 2099]
    seed::Int = 1234
    n_dhw_scenarios::Int = 50
    n_cyclone_scenarios::Int = 50
    reefs::ReefConfig = ReefConfig()
    connectivity::ConnectivityConfig = ConnectivityConfig()
    cover::CoverConfig = CoverConfig()
    dhw::DHWConfig = DHWConfig()
    cyclones::CycloneConfig = CycloneConfig()
end

"""Climate labels (sorted) for which a DHW file is written."""
climate_labels(cfg::DomainConfig)::Vector{Int} = sort!(collect(keys(cfg.dhw.climate)))

"""Simulation years covered by the domain."""
years(cfg::DomainConfig)::UnitRange{Int} = cfg.timeframe[1]:cfg.timeframe[2]

"""
    climate_params(dhw::DHWConfig, label::Int)::NamedTuple

Effective DHW parameters of one climate level (the `[dhw]` values with that level's
overrides applied).
"""
function climate_params(dc::DHWConfig, label::Int)::NamedTuple
    overrides = dc.climate[label]
    get_param(key) = haskey(overrides, key) ? overrides[key] : getfield(dc, Symbol(key))
    return (
        baseline=Float64(get_param("baseline")),
        trend_per_decade=Float64(get_param("trend_per_decade")),
        interannual_sd=Float64(get_param("interannual_sd")),
        ar1=Float64(get_param("ar1")),
        event_prob=Float64(get_param("event_prob")),
        event_magnitude=Float64.(collect(get_param("event_magnitude"))),
        event_years=Int.(collect(get_param("event_years"))),
        local_sd=Float64(get_param("local_sd"))
    )
end

# --------------------------------------------------------------------------------------
# Per-reef attribute specifications
# --------------------------------------------------------------------------------------

"""
    attr_spec(x)

Normalise a per-reef attribute specification to a `Float64`, a `Vector{Float64}`,
`(min=, max=)` or `(from=, to=)`.
"""
function attr_spec(x)
    x isa Real && return Float64(x)
    x isa AbstractVector && all(v -> v isa Real, x) && return Float64.(collect(x))

    if x isa AbstractDict || x isa NamedTuple
        d = Dict(String(k) => v for (k, v) in pairs(x))
        if keys(d) == Set(["min", "max"])
            return (min=Float64(d["min"]), max=Float64(d["max"]))
        elseif keys(d) == Set(["from", "to"])
            return (from=Float64(d["from"]), to=Float64(d["to"]))
        end
    end

    throw(ArgumentError(
        "Invalid per-reef attribute $(repr(x)): use a number, a list with one value per " *
        "reef, {min = a, max = b} or {from = a, to = b}."
    ))
end

"""
    resolve_attr(spec, n::Int, rng)::Vector{Float64}

Expand a per-reef attribute specification to one value per reef. Random ranges draw
`n` values sequentially from `rng`, so adding reefs does not change existing reefs.
"""
function resolve_attr(spec, n::Int, rng)::Vector{Float64}
    s = attr_spec(spec)
    if s isa Float64
        return fill(s, n)
    elseif s isa Vector{Float64}
        length(s) == n || throw(ArgumentError(
            "Per-reef list has $(length(s)) values but there are $n reefs."
        ))
        return copy(s)
    elseif haskey(s, :min)
        return s.min .+ (s.max - s.min) .* rand(rng, n)
    else
        return n == 1 ? [(s.from + s.to) / 2] : collect(range(s.from, s.to; length=n))
    end
end

"""Bounds of a per-reef attribute specification (for validation and UI defaults)."""
function attr_bounds(spec)::Tuple{Float64,Float64}
    s = attr_spec(spec)
    s isa Float64 && return (s, s)
    s isa Vector{Float64} && return isempty(s) ? (NaN, NaN) : extrema(s)
    return haskey(s, :min) ? (s.min, s.max) : extrema((s.from, s.to))
end

"""Expand a per-reef categorical value (single value or one per reef)."""
function resolve_categorical(spec, n::Int, T::Type)::Vector{T}
    if spec isa AbstractVector
        length(spec) == n || throw(ArgumentError(
            "Per-reef list has $(length(spec)) values but there are $n reefs."
        ))
        return T.(spec)
    end
    return fill(T(spec), n)
end

# --------------------------------------------------------------------------------------
# Random number streams
# --------------------------------------------------------------------------------------

# FNV-1a: a fixed string hash (Base.hash may change between Julia versions).
function _fnv1a(s::AbstractString)::UInt64
    h = 0xcbf29ce484222325
    for b in codeunits(s)
        h = (h ⊻ b) * 0x00000100000001b3
    end
    return h
end

"""
    component_rng(seed, tag)

Independent, reproducible random stream for one component of the domain (`"area_m2"`,
`"dhw"`, `"cyclones"`, ...). Separate streams mean that changing one setting (say the
connectivity density) does not reshuffle unrelated draws (say the DHW series).
StableRNGs keeps streams identical across Julia versions.
"""
component_rng(seed::Integer, tag::AbstractString) = StableRNG(_fnv1a(string(seed, ":", tag)))

# --------------------------------------------------------------------------------------
# TOML conversion
# --------------------------------------------------------------------------------------

_toml_value(x::NamedTuple) = Dict{String,Any}(String(k) => v for (k, v) in pairs(x))
_toml_value(x::AbstractVector) = collect(x)
_toml_value(x::Symbol) = String(x)
_toml_value(x) = x

function _section_dict(x)::Dict{String,Any}
    return Dict{String,Any}(
        String(f) => _toml_value(getfield(x, f)) for f in fieldnames(typeof(x))
    )
end

"""
    config_to_dict(cfg::DomainConfig)::Dict{String,Any}

TOML-ready dictionary of a configuration (inverse of [`load_config`](@ref)).
"""
function config_to_dict(cfg::DomainConfig)::Dict{String,Any}
    dhw = _section_dict(cfg.dhw)
    dhw["climate"] = Dict{String,Any}(
        string(label) => Dict{String,Any}(k => _toml_value(v) for (k, v) in overrides)
        for (label, overrides) in cfg.dhw.climate
    )

    return Dict{String,Any}(
        "domain" => Dict{String,Any}(
            "name" => cfg.name,
            "timeframe" => collect(cfg.timeframe),
            "seed" => cfg.seed,
            "n_dhw_scenarios" => cfg.n_dhw_scenarios,
            "n_cyclone_scenarios" => cfg.n_cyclone_scenarios
        ),
        "reefs" => _section_dict(cfg.reefs),
        "connectivity" => _section_dict(cfg.connectivity),
        "cover" => _section_dict(cfg.cover),
        "dhw" => dhw,
        "cyclones" => _section_dict(cfg.cyclones)
    )
end

"""
    save_config(path::AbstractString, cfg::DomainConfig)

Write a configuration as TOML (reload it with [`load_config`](@ref)).
"""
function save_config(path::AbstractString, cfg::DomainConfig)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "# SyntheticReefDomain.jl configuration")
        println(io, "# Regenerate with: generate_domain(\"$(basename(path))\")")
        println(io)
        TOML.print(io, config_to_dict(cfg); sorted=true)
    end
    return path
end

function _check_keys(d::AbstractDict, allowed, context::AbstractString)
    unknown = setdiff(String.(collect(keys(d))), allowed)
    isempty(unknown) || throw(ArgumentError(
        "Unknown key(s) in $context: $(join(unknown, ", ")). Allowed: $(join(allowed, ", "))."
    ))
    return nothing
end

# Convert one TOML value to the declared field type of a config struct.
function _field_value(T::Type, v)
    T === Float64 && return Float64(v)
    T === Int && return Int(v)
    T === Bool && return Bool(v)
    T === String && return String(v)
    T === Vector{Float64} && return Float64[Float64(x) for x in v]
    T === Vector{Int} && return Int[Int(x) for x in v]
    # `Any` fields: per-reef attribute or categorical specification
    v isa AbstractDict && return attr_spec(v)
    v isa AbstractVector && return identity.(collect(v))  # narrow Vector{Any}
    return v
end

function _from_dict(::Type{T}, d::AbstractDict, context::AbstractString) where {T}
    names = String.(fieldnames(T))
    _check_keys(d, names, context)
    kwargs = Dict{Symbol,Any}()
    for (k, v) in d
        f = Symbol(k)
        kwargs[f] = _field_value(fieldtype(T, f), v)
    end
    return T(; kwargs...)
end

function _parse_climate(d::AbstractDict)::Dict{Int,Dict{String,Any}}
    climate = Dict{Int,Dict{String,Any}}()
    for (k, overrides) in d
        key = string(k)
        occursin(r"^[1-9][0-9]*$", key) || throw(ArgumentError(
            "Climate label '$key' must be a positive integer without leading zeros " *
            "(ADRIA stores the RCP label as a number), e.g. [dhw.climate.45]."
        ))
        overrides isa AbstractDict || throw(ArgumentError(
            "[dhw.climate.$key] must be a table of parameter overrides."
        ))
        _check_keys(overrides, collect(CLIMATE_KEYS), "[dhw.climate.$key]")
        climate[parse(Int, key)] = Dict{String,Any}(
            String(p) => _field_value(fieldtype(DHWConfig, Symbol(p)), v)
            for (p, v) in overrides
        )
    end
    return climate
end

"""
    load_config(path::AbstractString)::DomainConfig
    load_config(d::AbstractDict)::DomainConfig

Read a configuration from a TOML file (or a parsed TOML dictionary). Missing keys take
their default values; unknown keys raise an error so typos are not silently ignored.
"""
load_config(path::AbstractString)::DomainConfig = load_config(TOML.parsefile(path))
function load_config(d::AbstractDict)::DomainConfig
    sections = ["domain", "reefs", "connectivity", "cover", "dhw", "cyclones"]
    _check_keys(d, sections, "the top level")

    dom = get(d, "domain", Dict{String,Any}())
    domain_keys = ["name", "timeframe", "seed", "n_dhw_scenarios", "n_cyclone_scenarios"]
    _check_keys(dom, domain_keys, "[domain]")
    kwargs = Dict{Symbol,Any}(
        Symbol(k) => _field_value(fieldtype(DomainConfig, Symbol(k)), v) for (k, v) in dom
    )

    kwargs[:reefs] = _from_dict(ReefConfig, get(d, "reefs", Dict()), "[reefs]")
    kwargs[:connectivity] = _from_dict(
        ConnectivityConfig, get(d, "connectivity", Dict()), "[connectivity]"
    )
    kwargs[:cover] = _from_dict(CoverConfig, get(d, "cover", Dict()), "[cover]")
    kwargs[:cyclones] = _from_dict(CycloneConfig, get(d, "cyclones", Dict()), "[cyclones]")

    dhw = Dict{String,Any}(String(k) => v for (k, v) in get(d, "dhw", Dict()))
    climate = pop!(dhw, "climate", nothing)
    dhw_cfg = _from_dict(DHWConfig, dhw, "[dhw]")
    if !isnothing(climate)
        dhw_cfg = DHWConfig(;
            (f => getfield(dhw_cfg, f) for f in fieldnames(DHWConfig))...,
            climate=_parse_climate(climate)
        )
    end
    kwargs[:dhw] = dhw_cfg

    return DomainConfig(; kwargs...)
end

# --------------------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------------------

function _check_attr!(errs, label, spec, n; lo=-Inf, hi=Inf, lo_open=false)
    s = try
        attr_spec(spec)
    catch err
        push!(errs, "$label: $(sprint(showerror, err))")
        return nothing
    end
    if s isa Vector{Float64} && length(s) != n
        push!(errs, "$label: per-reef list has $(length(s)) values but n = $n.")
    end
    if s isa NamedTuple && haskey(s, :min) && s.min > s.max
        push!(errs, "$label: min ($(s.min)) > max ($(s.max)).")
    end
    vmin, vmax = attr_bounds(s)
    below = lo_open ? vmin <= lo : vmin < lo
    if below || vmax > hi
        bound = lo_open ? "($lo, $hi]" : "[$lo, $hi]"
        push!(errs, "$label: values must lie in $bound, got $vmin to $vmax.")
    end
    return nothing
end

function _check_range!(errs, label, x, lo, hi)
    (lo <= x <= hi) || push!(errs, "$label must lie in [$lo, $hi], got $x.")
    return nothing
end

function _check_years!(errs, label, yrs, timeframe)
    outside = filter(y -> !(timeframe[1] <= y <= timeframe[end]), yrs)
    isempty(outside) || push!(errs, "$label outside the timeframe: $(join(outside, ", ")).")
    return nothing
end

function _check_interval!(errs, label, v, lo, hi)
    if length(v) != 2
        push!(errs, "$label must be [min, max], got $(repr(v)).")
    elseif v[1] > v[2] || v[1] < lo || v[2] > hi
        push!(errs, "$label must satisfy $lo <= min <= max <= $hi, got $(repr(v)).")
    end
    return nothing
end

"""
    validate(cfg::DomainConfig)::DomainConfig

Check a configuration and throw an `ArgumentError` listing every problem found.
Returns `cfg` unchanged when it is valid.
"""
function validate(cfg::DomainConfig)::DomainConfig
    errs = String[]

    # [domain]
    occursin(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$", cfg.name) || push!(errs,
        "name '$(cfg.name)' must use only letters, digits, '_', '-' or '.' " *
        "(it becomes the folder and geopackage name)."
    )
    tf = cfg.timeframe
    if length(tf) != 2 || tf[1] >= tf[2]
        push!(errs, "timeframe must be [start_year, end_year] with start < end.")
    elseif tf[2] - tf[1] + 1 < MIN_YEARS
        push!(errs, "timeframe must span at least $MIN_YEARS years, got $(tf[2] - tf[1] + 1).")
    end
    cfg.n_dhw_scenarios >= 1 || push!(errs, "n_dhw_scenarios must be >= 1.")
    cfg.n_cyclone_scenarios >= 1 || push!(errs, "n_cyclone_scenarios must be >= 1.")

    # [reefs]
    r = cfg.reefs
    n = r.n
    1 <= n <= MAX_REEFS || push!(errs, "reefs.n must be between 1 and $MAX_REEFS, got $n.")
    r.layout in LAYOUTS ||
        push!(errs, "reefs.layout must be one of $(join(LAYOUTS, ", ")), got '$(r.layout)'.")
    if length(r.center) != 2
        push!(errs, "reefs.center must be [lon, lat].")
    else
        _check_range!(errs, "reefs.center longitude", r.center[1], -180.0, 180.0)
        _check_range!(errs, "reefs.center latitude", r.center[2], -80.0, 80.0)
    end
    r.spacing_km > 0 || push!(errs, "reefs.spacing_km must be > 0.")
    r.cluster_spacing_km > 0 || push!(errs, "reefs.cluster_spacing_km must be > 0.")
    1 <= r.n_clusters <= max(n, 1) ||
        push!(errs, "reefs.n_clusters must be between 1 and n ($n), got $(r.n_clusters).")
    _check_attr!(errs, "reefs.area_m2", r.area_m2, n; lo=0.0, lo_open=true)
    _check_attr!(errs, "reefs.k_percent", r.k_percent, n; lo=0.0, hi=100.0, lo_open=true)
    _check_attr!(errs, "reefs.depth_med", r.depth_med, n; lo=0.0, lo_open=true)
    _check_attr!(errs, "reefs.ub_med", r.ub_med, n; lo=0.0)
    for (label, spec, T) in (("reefs.habitat", r.habitat, String),
        ("reefs.cb_calib_group", r.cb_calib_group, Int))
        try
            resolve_categorical(spec, n, T)
        catch err
            push!(errs, "$label: $(sprint(showerror, err))")
        end
    end
    groups = r.cb_calib_group isa AbstractVector ? r.cb_calib_group : [r.cb_calib_group]
    all(g -> g isa Integer && 1 <= g <= 12, groups) ||
        push!(errs, "reefs.cb_calib_group values must be integers in 1-12.")

    # [connectivity]
    c = cfg.connectivity
    c.topology in TOPOLOGIES || push!(errs,
        "connectivity.topology must be one of $(join(TOPOLOGIES, ", ")), got '$(c.topology)'."
    )
    _check_range!(errs, "connectivity.self_recruitment", c.self_recruitment, 0.0, 1.0)
    _check_range!(errs, "connectivity.strength", c.strength, 0.0, 1.0)
    _check_range!(errs, "connectivity.asymmetry", c.asymmetry, 0.0, 1.0)
    _check_range!(errs, "connectivity.density", c.density, 0.0, 1.0)
    c.reach >= 1 || push!(errs, "connectivity.reach must be >= 1.")
    c.decay_km > 0 || push!(errs, "connectivity.decay_km must be > 0.")

    # [cover]
    cv = cfg.cover
    cv.basis in COVER_BASES ||
        push!(errs, "cover.basis must be one of $(join(COVER_BASES, ", ")), got '$(cv.basis)'.")
    _check_attr!(errs, "cover.total_percent", cv.total_percent, n; lo=0.0, hi=100.0)
    if length(cv.composition) != N_GROUPS
        push!(errs, "cover.composition needs $N_GROUPS weights ($(join(FUNCTIONAL_GROUPS, ", "))).")
    elseif any(<(0), cv.composition) || sum(cv.composition) <= 0
        push!(errs, "cover.composition weights must be >= 0 with a positive sum.")
    end
    cv.jitter >= 0 || push!(errs, "cover.jitter must be >= 0.")

    # [dhw]
    d = cfg.dhw
    isempty(d.climate) && push!(errs, "dhw.climate needs at least one climate level.")
    for label in sort(collect(keys(d.climate)))
        ctx = "dhw.climate.$label"
        label >= 1 || push!(errs, "$ctx: climate labels must be positive integers.")
        unknown = setdiff(keys(d.climate[label]), CLIMATE_KEYS)
        if !isempty(unknown)
            push!(errs, "$ctx: unknown override(s) $(join(unknown, ", ")).")
            continue
        end
        p = climate_params(d, label)
        p.interannual_sd >= 0 || push!(errs, "$ctx: interannual_sd must be >= 0.")
        p.local_sd >= 0 || push!(errs, "$ctx: local_sd must be >= 0.")
        (0 <= p.ar1 < 1) || push!(errs, "$ctx: ar1 must lie in [0, 1).")
        _check_range!(errs, "$ctx: event_prob", p.event_prob, 0.0, 1.0)
        _check_interval!(errs, "$ctx: event_magnitude", p.event_magnitude, 0.0, Inf)
        length(tf) == 2 && _check_years!(errs, "$ctx: event_years", p.event_years, tf)
    end
    _check_attr!(errs, "dhw.reef_exposure", d.reef_exposure, n; lo=0.0)
    _check_attr!(errs, "dhw.reef_offset", d.reef_offset, n)

    # [cyclones]
    cy = cfg.cyclones
    _check_range!(errs, "cyclones.annual_prob", cy.annual_prob, 0.0, 1.0)
    _check_interval!(errs, "cyclones.intensity", cy.intensity, 0.0, 1.0)
    if length(cy.group_vulnerability) != N_GROUPS
        push!(errs, "cyclones.group_vulnerability needs $N_GROUPS values.")
    elseif any(v -> !(0 <= v <= 1), cy.group_vulnerability)
        push!(errs, "cyclones.group_vulnerability values must lie in [0, 1].")
    end
    length(tf) == 2 && _check_years!(errs, "cyclones.event_years", cy.event_years, tf)

    isempty(errs) || throw(ArgumentError(
        "Invalid synthetic domain configuration:\n  - " * join(errs, "\n  - ")
    ))
    return cfg
end
