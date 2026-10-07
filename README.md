# SyntheticReefDomain.jl

Generates synthetic reef systems as ADRIA domain data packages. You choose the number of
reefs (1-20), how they are connected, their initial coral cover and the disturbances they
experience (heat stress and cyclones). The output folder loads with
`ADRIA.load_domain(path, 45)` like any real domain, so it also runs through every engine that
reads an ADRIA `Domain` (CoralBlox, CoralPop, Kora, ReefMod, CORSET, and C~scape via
CscapeInterface). Use it for controlled experiments that real domains cannot give you: a
handful of reefs, known source-sink structure, refuges vs. hotspots, cyclone-heavy regimes, ...

There are two ways in:

- **Script / package**: describe the domain in a TOML file (or with keyword arguments) and
  call `generate_domain`.
- **Notebook**: a Pluto notebook with sliders and live previews (map, connectivity, cover,
  DHW and cyclone time series) that writes the domain or saves its TOML config.

## Status

- Output verified against ADRIA 0.17.2 (`ADRIA.jl` branch `CscapeIntegration`).
  - Domains with 1, 2, 10 and 20 reefs and every connectivity topology load, and their
    connectivity, DHW, waves, cyclone mortality and per-group initial cover come back
    identical.
  - Counterfactual scenarios run for 1-20 reefs and for several climate levels in one
    `run_scenarios` call. Intervention-sampled runs (`ADRIA.sample`) were checked with 1, 2
    and 10 reefs.
- C~scape translation (`CscapeInterface.domain_to_cscape_inputs`) picks up `ub_med` from the
  geopackage (`sandbox/02_cscape_translate.jl`).
- Not covered: waves as a disturbance (ADRIA's wave mortality is switched off), 5-D MCB DHW
  data, and multiple sites per reef (each reef is one ADRIA location).

## Setup

Julia 1.11. The test, sandbox and notebook environments find ADRIA, ADRIAIndicators and this
package through `[sources]` paths, so keep the sibling layout:

```
GitHub/
  ADRIA.jl/            (ADRIA package in ADRIA.jl/ADRIA)
  ADRIAIndicators.jl/
  SyntheticReefDomain.jl/
```

```julia
julia --project=. -e "using Pkg; Pkg.instantiate()"            # the package itself
julia --project=test -e "using Pkg; Pkg.instantiate()"         # tests (pulls in ADRIA)
julia --project=notebooks -e "using Pkg; Pkg.instantiate()"    # Pluto notebook
julia --project=sandbox -e "using Pkg; Pkg.instantiate()"      # example scripts
```

## Quick start

```julia
using SyntheticReefDomain

# From a config file (writes domains/SyntheticReefs)
path = generate_domain("configs/default.toml")

# Or in code
cfg = DomainConfig(;
    name="Chain5",
    reefs=ReefConfig(; n=5, layout="line"),
    connectivity=ConnectivityConfig(; topology="chain", asymmetry=1.0),
    cover=CoverConfig(; total_percent=(min=10.0, max=30.0)),
    dhw=DHWConfig(; reef_exposure=(from=0.6, to=1.4)),
    cyclones=CycloneConfig(; annual_prob=0.2)
)
dom = build_domain(cfg)          # in memory, nothing written; `display(dom)` summarises it
write_domain(dom, "domains/Chain5")
save_config("configs/Chain5.toml", cfg)

# In ADRIA
using ADRIA
adria_dom = ADRIA.load_domain("domains/Chain5", 45)
```

`configs/default.toml` documents every option with its default value. The examples in
`configs/examples/` show a downstream-only chain with three climate levels, isolated reefs with
explicit refuges, and a cyclone-heavy clustered system.

## What you can set

| Section | Settings |
| --- | --- |
| `[domain]` | name, timeframe (at least 26 years), seed, number of DHW and cyclone scenarios |
| `[reefs]` | number of reefs (1-20); layout (`line`, `grid`, `ring`, `random`, `clustered`); spacing; bearing; clusters; area, habitable fraction `k`, depth, `ub_med`, habitat, CoralBlox calibration group |
| `[connectivity]` | topology (`isolated`, `chain`, `ring`, `all_to_all`, `distance_decay`, `random`); self-recruitment; strength; asymmetry and reach (chain/ring); decay distance; link density |
| `[cover]` | total cover per reef, as % of habitable or of whole-reef area; functional-group composition; per-reef variation of the mix |
| `[dhw]` | baseline; warming trend; year-to-year variability (SD, autocorrelation); marine heatwaves (probability, size, forced years); per-reef exposure and offset; local noise; one or more climate levels |
| `[cyclones]` | on/off, annual probability, forced years, intensity range, mortality per functional group |

Per-reef attributes (area, k, depth, `ub_med`, total cover, DHW exposure and offset) take any
of four forms:

| Form | TOML | Julia |
| --- | --- | --- |
| same value on every reef | `k_percent = 50.0` | `50.0` |
| one value per reef | `k_percent = [40.0, 55.0, 60.0]` | `[40.0, 55.0, 60.0]` |
| uniform random range | `k_percent = {min = 45.0, max = 63.0}` | `(min=45.0, max=63.0)` |
| gradient, first to last reef | `reef_exposure = {from = 0.6, to = 1.4}` | `(from=0.6, to=1.4)` |

Reef order follows the layout: west to east along a line with bearing 90, serpentine rows
for a grid, around the ring, cluster by cluster. Gradients and chain connectivity follow
that order.

### Heat stress

Each scenario has a regional DHW signal: baseline, plus a linear trend, plus an AR(1)
anomaly, plus random or forced marine heatwaves. Each reef scales and shifts that signal and
adds its own noise:

```
DHW[t, reef, s] = max(0, exposure[reef] * R[t, s] + offset[reef] + local_sd * noise)
```

So reefs share bad years but differ in severity, as in ADRIA's Test_domain (between-reef
correlation about 0.4-0.75). `local_sd = 0` with exposure 1 gives identical DHW everywhere.
The defaults give about 5 DHW in the first decade and 12 in the last, close to Test_domain.

### Climate levels ("RCPs")

`[dhw.climate.<label>]` tables each write a `DHWs/dhwRCP<label>.nc`. A table can override
the trend, baseline, variability, heatwave settings or local noise for that level. ADRIA
treats the label purely as a number: `run_scenarios` parses it with `parse(Int, ...)`, and
`load_domain(path, 45)` picks the file. Any positive integer works (1, 2, 3 or 26, 45, 85), and
it says nothing about real emission pathways. The default is a single level, `45`, so
existing scripts that load domains with RCP 45 work unchanged. All levels share the same
random draws: scenario 7 under label 1 and under label 2 differ only by their parameters.

### Cyclones

Cyclones are system-wide. In a cyclone year every reef loses the same fraction of each
functional group: `intensity * group_vulnerability[g]`. By default tabular Acropora are the
most vulnerable and large massives the least.

### Waves

There is no wave disturbance (ADRIA's wave mortality is commented out). The wave file holds
each reef's `ub_med` as a constant series, because:

- ADRIA's site selection uses waves, and indexes them with the DHW scenario, so the file needs
  as many scenarios as the DHW file;
- C~scape reads `ub_med` from the geopackage.

## Notebook

```julia
julia --project=notebooks -e "using Pluto; Pluto.run(notebook=\"notebooks/reef_builder.jl\")"
```

1. Pick a starting config (`configs/*.toml` and `configs/examples/*.toml`).
2. Adjust the controls. The domain is rebuilt in memory on every change and previewed:
   - a reef map with connectivity arrows;
   - a connectivity matrix;
   - cover by functional group;
   - DHW for one scenario, plus the spread across scenarios and climate levels;
   - cyclone years and mortality per group.
3. Press **Write domain** (default `domains/<name>`) and/or **Save config** (default
   `configs/<name>.toml`). Neither overwrites existing files unless *replace* is ticked.

Per-reef lists from the starting file are kept when their length matches the number of
reefs. Untick the checkbox at the top to use the range controls instead.

## Running in ADRIA

```julia
julia --project=sandbox -t 4 sandbox/01_generate_and_run.jl configs/examples/chain_source_sink.toml
```

This generates the domain, runs 16 counterfactual scenarios for every climate level, and
prints mean cover per decade. `sandbox/02_cscape_translate.jl` checks the C~scape
translation; it runs in CscapeInterface's own environment.

Remember to call `ADRIA.setup()` before `run_scenarios`, as in the scripts. Without it ADRIA
fails with `KeyError: ADRIA_THRESHOLD`.

## Domain format

The output mirrors `ADRIA.jl/ADRIA/test/data/Test_domain` (data package version 0.8.0):

| Path | Content |
| --- | --- |
| `datapackage.json` | `version` 0.8.0, `simulation_metadata.timeframe`, resource list |
| `spatial/<Name>.gpkg` | EPSG:4326 reef polygons. Columns: `reef_siteid`/`UNIQUE_ID` (`reef_01`...), `cluster_id`, `k` (%), `area` (m²), `depth_med`, `ub_med`, `CB_CALIB_GROUPS`, plus informational `habitat`, `zone_type`, `dhw_exposure`, `dhw_offset`, `lon`, `lat`, `x_km`, `y_km` |
| `spatial/coral_cover.nc` | `layer(species=5, locations)`, fraction of habitable area; ADRIA splits it into 7 size classes |
| `connectivity/connectivity_<Y0>_synthetic.csv` | source (rows) x sink (columns) transfer probabilities |
| `DHWs/dhwRCP<label>.nc` | `dhw(timesteps, sites, member)` and a `reef_siteid` variable |
| `waves/wave_RCP<label>.nc` | `Ub(timesteps, sites, member)` |
| `cyclones/cyclone_mortality.nc` | `layer(timesteps, locations, species=5, scenarios)` |
| `generator_config.toml`, `README.md` | the exact config used (regenerate with `generate_domain`) and a summary |

## Things to know

- **ADRIA names a domain after its folder.** The geopackage must be called
  `<folder name>.gpkg`; `write_domain` handles this.
- **Connectivity in ADRIA's own model is relative.**
  - ADRIA divides each sink's column by its sum, so a reef's larval supply is a weighted mix of
    its sources' output. Strength changes *which* reefs supply a reef, not *how much* arrives.
  - The one way to cut a reef off entirely is no incoming links at all, i.e. `isolated` with
    `self_recruitment = 0`. In a quick test, the default domain's mean cover rose from 0.29 to
    about 0.5, but with no connectivity it fell to 0.01 by the 2090s.
  - Absolute values do matter for ADRIA's connectivity-based decision criteria and for engines
    that move larvae through the matrix (CoralPop, Kora, C~scape).
- **One reef works but has no neighbours.** ADRIA's `mean_to_neighbor` is NaN for a single
  location. Counterfactual and intervention runs still complete.
- **Timeframe of at least 26 years.**
  - ADRIA's intervention start-year factors go up to year 25.
  - ADRIA's NetCDF loader finds the location dimension by length, so the number of years must
    differ from the number of reefs. With at most 20 reefs, 26 years guarantees that.
- **Identical `(ub_med, depth_med)` pairs** make C~scape merge reefs onto one demography
  point. Exact duplicates are nudged apart by 1e-4 m/s unless `distinct_ub_depth = false`.
- **Reproducibility.** Every component draws from its own seeded StableRNG stream. So:
  - the same config gives the same domain across Julia versions (except cover `jitter`,
    whose Dirichlet draws follow the installed Distributions.jl);
  - changing one section (say connectivity density) does not reshuffle another (say the DHW
    series);
  - adding reefs leaves the existing reefs' attributes and noise unchanged.
- Reef polygons are irregular ellipses with exactly the requested area. At the default
  spacing they are small dots on the map, but ADRIA only uses their centroids.

## Tests

```julia
julia --project=test -t 4 test/runtests.jl
```

The tests cover:

- configuration parsing and validation;
- every layout and topology;
- cover conversion;
- disturbance statistics;
- an ADRIA load round-trip for 1-20 reefs;
- short ADRIA scenario runs.

Set `RUN_ADRIA_SCENARIOS=false` to skip the scenario runs, which take about 2-3 minutes.
