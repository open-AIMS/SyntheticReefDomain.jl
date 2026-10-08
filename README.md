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

A second, separate notebook exports a domain as an input project for the original R package
[cscape](https://github.com/open-AIMS/c-scape) (see *C~scape (R) export*).

## Getting started

**With Julia (Windows):**
1. Install Julia 1.11 or newer (https://julialang.org/install).
2. Get this repository.
3. Double-click a launcher:
   - `start_reef_builder.bat` for the domain builder;
   - `start_cscape_export.bat` for the C~scape export.

The first start installs and compiles the packages, which takes 10-15 minutes. Later starts
take under a minute. The notebook opens in your web browser. The first time, Pluto may show
*Safe preview*: click **Run notebook code**.

**Without Julia (Binder):** once this repository is public on GitHub, the notebooks run in
the browser through [mybinder.org](https://mybinder.org). It is free and needs no account or
install. Replace `OWNER/REPO` with the GitHub location:

- Domain builder:
  `https://mybinder.org/v2/gh/OWNER/REPO/HEAD?urlpath=pluto%2Fopen%3Fpath%3Dnotebooks%2Freef_builder.jl`
- C~scape export:
  `https://mybinder.org/v2/gh/OWNER/REPO/HEAD?urlpath=pluto%2Fopen%3Fpath%3Dnotebooks%2Fcscape_export.jl`

Binder limits:
- the first launch builds the image, which is slow; later launches are faster;
- about 1-2 GB of memory;
- idle sessions stop after about 10 minutes, and everything in a session is deleted when it
  ends. Use **Prepare download (.zip)** to keep your domain or C~scape project;
- no R, so the C~scape export there contains the inputs plus a script to run in R later.

The Binder setup lives in `binder/` (following
[pluto-on-binder](https://github.com/fonsp/pluto-on-binder)):
- `Project.toml` picks the Julia version;
- `requirements.txt` and `jupyter_server_config.py` serve Pluto through jupyter-server-proxy;
- `postBuild` installs and precompiles the notebook environment and opens both notebooks
  once (`warmup.jl`).

Before making the repository public, check the Moore-derived default ranges with their
owners (see *Provenance*).

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
- C~scape (R) export verified with cscape 1.0.0 (R 4.4.2). An exported project builds with
  R, passes `run_cscape(..., dry_run = TRUE)`, and runs a full replicate (75 years) with the
  dummy demography of CScape.jl's `00_Verify` project.
- The Binder pieces were checked locally (warm-up and the Pluto start command). The
  mybinder.org launch itself can only be tested once the repository is public.
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
using CairoMakie   # optional: adds the plots to each domain's README.md

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

### Domain report

Every written domain gets a `README.md` that documents it:

- a summary and the `ADRIA.load_domain` call;
- every setting, in the notebook's sections (Domain, Reefs, Connectivity, Initial cover,
  Heat stress, Cyclones), including per-climate-level overrides;
- a per-reef table of the resolved values (position, area, k, depth, `ub_med`, DHW exposure
  and offset, initial cover, mean DHW);
- the notebook's plots, saved as PNGs in `figures/`: reef map with connectivity,
  connectivity matrix, initial cover, DHW per climate level and cyclones;
- a list of the files in the package.

The plots come from an optional extension and are included when CairoMakie is loaded
(`using CairoMakie`; the notebook always loads it). Without it, or with
`write_domain(dom, path; figures=false)`, the README still lists every setting and says how
to add the plots. The same plots are available directly as `plot_map(dom)`,
`plot_connectivity(dom)`, `plot_cover(dom)`, `plot_dhw(dom, label, scenario)` and
`plot_cyclones(dom, scenario)`.

## What you can set

| Section | Settings |
| --- | --- |
| `[domain]` | name, timeframe (at least 26 years), seed, number of stochastic DHW and cyclone scenarios (scenario 1 is always disturbance-free; see below) |
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

### Scenario 1: no disturbance

Every domain has a disturbance-free control scenario. Scenario 1 of every DHW file (each
climate level) and of the cyclone file is all zeros, for every reef and year. The configured
stochastic scenarios follow as 2..n+1, so with `n_dhw_scenarios = 50` the DHW files hold 51
scenarios. The wave file matches that count; its scenario 1 holds the usual constant
`ub_med`, since waves are not a disturbance here.

In ADRIA:

- `dhw_scenario` / `cyclone_mortality_scenario` = 1 selects the control scenario. It is also
  these factors' nominal value, so runs that use nominal factor values are disturbance-free.
- Sampling (`ADRIA.sample`, `sample_cf`) draws uniformly over all scenarios, so about 1 in
  n+1 sampled runs is disturbance-free. To leave it out:

  ```julia
  n = size(dom.dhw_scens, 3)                      # stochastic + 1
  ADRIA.set_factor_bounds!(dom, :dhw_scenario, Tuple(2.0:Float64(n)))
  ADRIA.set_factor_bounds!(dom, :cyclone_mortality_scenario,
      Tuple(2.0:Float64(size(dom.cyclone_mortality_scens, 4))))
  ```

### Waves

There is no wave disturbance (ADRIA's wave mortality is commented out). The wave file holds
each reef's `ub_med` as a constant series, because:

- ADRIA's site selection uses waves, and indexes them with the DHW scenario, so the file needs
  as many scenarios as the DHW file;
- C~scape reads `ub_med` from the geopackage.

## Notebook

```
julia --project=notebooks notebooks/run.jl
```

(From the repo root. Or, inside Julia: `using Pluto; Pluto.run(notebook="notebooks/reef_builder.jl")`.)

1. Pick a starting config. You can do this in two ways:
   - **From the dropdown**, which lists `configs/*.toml`, `configs/examples/*.toml` and the
     `generator_config.toml` of every domain already in `domains/`, so you can continue from
     a domain you wrote earlier. The list is read when the notebook starts.
   - **Open any `.toml` file** with the file button, which opens the Windows file dialog.
     A file opened this way wins over the dropdown. Choosing from the dropdown again clears
     it. If the file is not a valid config, a warning says why and the dropdown choice is
     used instead.
2. Adjust the controls. The domain is rebuilt in memory on every change and previewed:
   - a reef map with connectivity arrows;
   - a connectivity matrix;
   - cover by functional group;
   - DHW for one scenario, plus the spread across scenarios and climate levels;
   - cyclone years and mortality per group.
3. Press **Write domain** (default `domains/<name>`) and/or **Save config** (default
   `configs/<name>.toml`). Neither overwrites existing files unless *replace* is ticked.
   The written domain's `README.md` records every setting and shows the same plots (see
   *Domain report*). **Prepare download (.zip)** gives the same domain as a zip, which is
   useful on Binder.

Per-reef lists from the starting file are kept when their length matches the number of
reefs. Untick the checkbox at the top to use the range controls instead.

## C~scape (R) export

`notebooks/cscape_export.jl` (or `start_cscape_export.bat`, or
`julia --project=notebooks notebooks/run_cscape.jl`) turns a synthetic domain into an input
project for the R package [cscape](https://github.com/open-AIMS/c-scape). It is kept
separate from the domain builder.

1. Choose the domain: one already in `domains/`, a config, or any `.toml` file. The domain
   is rebuilt from its settings, so the C~scape inputs match the ADRIA domain exactly.
2. Set the C~scape options:
   - Region: a `Cluster` level of the demography models, e.g. `Offshore_Central`.
   - Demography file names.
   - Spawner connectivity scale.
   - Brooder dispersal.

   The options of an earlier export of the same domain are loaded as the starting point.
3. Check the previews: the site table, spawner and brooder connectivity, cyclone
   categories per replicate, and initial cover for the six C~scape types. Warnings flag
   reefs outside the demography models' depth or `ub_med` range.
4. **Write C~scape project** (default `<domain>/cscape/`) or **Prepare download (.zip)**.

From Julia:

```julia
write_cscape_project(dom, "domains/Chain5/cscape"; options=CscapeOptions(; region="Offshore_Central"))
```

Julia cannot write `.rds`, so the export writes CSV inputs (`inputs/`) plus
`build_cscape_project.R`. The script turns them into `ScenarioID.xlsx` and
`data/*.rds` using cscape's own helpers. It runs automatically when R and cscape are
installed; otherwise run `Rscript build_cscape_project.R` later.

**Demography files are not part of the export.** Copy the six `demog_inputs_<type>.rds`
files into the project's `data/` folder. Then run in R:
`cscape::run_cscape_batch(45, "<project folder>")`. Scenario IDs are the climate labels.

| C~scape input | From the synthetic domain |
| --- | --- |
| site table (`reef_sites_<name>.rds`) | `reef_siteid`, `k` (%), `area` (m²), `depth_med`, `ub_med`, plus `long`/`lat`, `habitat`, `cluster_id` |
| connectivity `[destination, source, type]` | Spawning types: the domain matrix transposed (x scale). Brooder: 90% of `brood_site_value` stays home, 10% goes to reefs within the radius (as in CscapeInterface) |
| temporal file, one per climate label | Replicate r = DHW scenario r, so **replicate 1 is the disturbance-free control**. `DHW_site` = per-reef DHW. Cyclone = category `ceil(5 x intensity)` (1-5), the same on every reef. No COTS, no external eggs, `temp_growth` 1 |
| `initialisationData.csv` | % of site area. Corymbose non-Acropora is split evenly between `corym_non_acro` and `brooder` |
| `ScenarioID.xlsx` | one row per climate label, `num_repeats` = all replicates, `Region` and `Growth_Surv_file` from the options |

Differences from the ADRIA run to keep in mind:
- **Cyclones:** C~scape computes cyclone mortality from the category (wind speed), so the
  domain's per-group cyclone vulnerabilities are not used.
- **Connectivity:** C~scape uses connectivity as absolute transfer probabilities, while
  ADRIA's own model only uses proportions. Use the scale option if the domain's values are
  too high for C~scape.

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
| `DHWs/dhwRCP<label>.nc` | `dhw(timesteps, sites, member)` and a `reef_siteid` variable; member 1 is all zeros |
| `waves/wave_RCP<label>.nc` | `Ub(timesteps, sites, member)` |
| `cyclones/cyclone_mortality.nc` | `layer(timesteps, locations, species=5, scenarios)`; scenario 1 is all zeros |
| `generator_config.toml` | the exact config used (regenerate with `generate_domain`) |
| `README.md`, `figures/*.png` | report of all settings, a per-reef table and the plots (see *Domain report*) |

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

## Provenance

All code in this repository was written for it (with Claude Code, October 2026). No code was
copied from other projects, papers or websites. Other projects were used in only two ways:
their code and data were read to work out the format ADRIA expects and to pick default
values, and open-source libraries are used as ordinary dependencies.

**None of the ecological default values are peer reviewed or calibrated.** They are
placeholders for building plausible test systems. Replace them with literature values, and
cite those, before using results in a publication.

### Format

| Source | Used for | Licence | Peer reviewed |
| --- | --- | --- | --- |
| [ADRIA.jl](https://github.com/open-AIMS/ADRIA.jl) (AIMS, DOI [10.25845/E9VT-4G35](https://doi.org/10.25845/E9VT-4G35)): loaders in `Domain.jl`, `io/inputs.jl`, `io/initial_coral_cover.jl`, `ecosystem/connectivity.jl`, `scenario.jl` | File names, NetCDF layouts, required columns, ID handling, RCP label parsing, connectivity normalisation | MIT | No (software archive) |
| ADRIA `test/data/Test_domain` | Layout of every output file | MIT | No (test data) |
| [CscapeInterface.jl](https://github.com/open-AIMS/CscapeAccess) `domain_translator.jl` | The `ub_med` column C~scape reads | MIT | No |
| CScape.jl `KNOWN_QUIRKS.md` (documentation only, no code) | Distinct `(ub_med, depth_med)` pairs | GPL | No |
| [cscape](https://github.com/open-AIMS/c-scape) R package: `run_cscape.R`, `load_data.R`, `site_loop.R`, `ipm_pred.R`, `simulate_reef.R` (read, not copied) | C~scape project layout, file structures, cyclone categories, scenario-table columns. The generated `build_cscape_project.R` *calls* cscape's exported functions | GPL-3 | No |
| CscapeInterface.jl `domain_translator.jl` | Conventions for the C~scape export: transposed connectivity, brooder matrix, cover split, units | MIT | No |
| [pluto-on-binder](https://github.com/fonsp/pluto-on-binder) | Approach for `binder/` (jupyter-server-proxy, Pluto settings) | MIT | No |
| [open-AIMS/ADRIA-data-specs](https://github.com/open-AIMS/ADRIA-data-specs) | Not consulted. This is the authoritative spec to check against. | -- | -- |

### Default values

| Setting | Default | Source | Status |
| --- | --- | --- | --- |
| `area_m2`, `k_percent`, `depth_med`, `ub_med` | 60,000-93,000 m², 45-63 %, 2-15 m, 0.23-0.55 m/s | Interquartile ranges of the Moore_2025-08-22_v080 ADRIA domain (AIMS / RRAP; spatial data from A. Cresswell) | No licence stated, unpublished: confirm with the data owners before a public release |
| Habitat names | Crest, Outer Flat, Sheltered Slope, Slope | Moore domain categories | As above |
| `cover.composition` | 0.26, 0.25, 0.07, 0.19, 0.23 | Mean cover per functional group in ADRIA's Test_domain | Test data, not observations |
| `cover.total_percent` | 10-40 % of habitable area | Placeholder (Test_domain per-location totals span about 2-45 %) | Not calibrated |
| `[dhw]` parameters | baseline 3.5, trend 1.2 DHW/decade, SD 2, AR(1) 0.3, heatwave probability 0.15, size 2-6 DHW, exposure 0.7-1.3, local SD 2 | Tuned so the default domain resembles Test_domain DHW (about 5 DHW in the first decade, 13 in the last) | Not calibrated to observations or climate projections |
| `[cyclones]` parameters | annual probability 0.1, intensity 0.2-1.0, vulnerability 0.6 / 0.5 / 0.4 / 0.2 / 0.1 | Placeholders. The ordering (Acropora most, massives least susceptible) follows general understanding; the numbers come from no source | Not calibrated |
| `[connectivity]` parameters | distance decay over 10 km, self-recruitment 0.2, strength 0.05 | Illustrative placeholders | Not calibrated |
| `reefs.center` | 146.0° E, 18.0° S | Arbitrary location in the central Great Barrier Reef region | -- |

### Methods

The methods are standard, written from general knowledge without a specific source:

- great-circle (haversine) distances;
- an equirectangular km to lon/lat approximation;
- shoelace polygon areas;
- an AR(1) anomaly process;
- Dirichlet sampling of cover composition;
- an FNV-1a hash that derives the random-number seed for each part of the domain.

The model structures are the author's own simple choices, not taken from the literature:

- reef DHW = exposure x regional signal + offset + noise;
- cyclone mortality = intensity x group vulnerability;
- the connectivity topologies.

### Dependencies

All used as libraries; no code was copied from them.

| Package | Licence |
| --- | --- |
| NetCDF.jl, GeoDataFrames.jl, ArchGDAL.jl, GeoFormatTypes.jl, DataFrames.jl, Distributions.jl, StableRNGs.jl, JSON.jl, ZipArchives.jl | MIT |
| R packages called by the generated C~scape build script (cscape, sf, writexl); not shipped here | GPL-3 (cscape), GPL-2 or MIT (sf), BSD-2 (writexl) |
| Pluto.jl (notebook only) | MIT |
| CairoMakie.jl (notebook, and optional plots in the domain README via a package extension) | MIT |
| PlutoUI.jl (notebook only) | Unlicense (public domain) |
| ADRIA.jl, ADRIAIndicators.jl (tests and sandbox only) | MIT |

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
