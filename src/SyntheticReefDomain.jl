"""
    SyntheticReefDomain

Generate synthetic ADRIA domain data packages: a system of 1-20 reefs with a configurable
layout, connectivity, initial coral cover and disturbances (heat stress and cyclones).

The output folder is loadable with `ADRIA.load_domain(path, <climate label>)` and therefore
by every engine that reads an ADRIA `Domain` (ADRIA/CoralBlox, CoralPop, Kora, ReefMod,
CORSET, C~scape via CscapeInterface).

| ADRIA domain file                    | Built by              | Written by             |
| ------------------------------------ | --------------------- | ---------------------- |
| `spatial/<Name>.gpkg`                | `build_reefs`         | `write_gpkg`           |
| `connectivity/connectivity_<Y0>_*.csv` | `build_connectivity` | `write_connectivity`   |
| `spatial/coral_cover.nc`             | `build_cover`         | `write_cover_nc`       |
| `DHWs/dhwRCP<label>.nc`              | `build_dhw`           | `write_env_nc`         |
| `waves/wave_RCP<label>.nc`           | `build_waves`         | `write_env_nc`         |
| `cyclones/cyclone_mortality.nc`      | `build_cyclones`      | `write_cyclone_nc`     |
| `datapackage.json`                   | --                    | `write_datapackage`    |

The format mirrors `ADRIA.jl/ADRIA/test/data/Test_domain` (data package version 0.8.0);
see the README for the full contract and the reasons behind each constraint.
"""
module SyntheticReefDomain

using Dates, Printf, Random, Statistics, TOML

using DataFrames
using Distributions: Dirichlet
using StableRNGs: StableRNG

import ArchGDAL as AG
import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import JSON
import NetCDF
import ZipArchives

# Configuration: structs, TOML (de)serialisation and validation
include("config.jl")

# Builders: pure functions of the config (seeded, no file I/O)
include("layout.jl")
include("connectivity.jl")
include("cover.jl")
include("disturbances.jl")

# Assembly and output
include("generate.jl")
include("writers.jl")
include("report.jl")

# Export to the original R package cscape (separate notebook: notebooks/cscape_export.jl)
include("cscape.jl")

# Plots: methods are added by the CairoMakie extension (ext/), loaded with `using CairoMakie`
"""
    plot_map(dom::SyntheticDomain)

Reef positions (coloured by initial cover, sized by area, outlined by cluster) with
connectivity arrows. Requires `using CairoMakie`.
"""
function plot_map end

"""
    plot_connectivity(dom::SyntheticDomain)

Source x sink connectivity matrix as a heatmap. Requires `using CairoMakie`.
"""
function plot_connectivity end

"""
    plot_cover(dom::SyntheticDomain)

Initial cover per reef, stacked by functional group. Requires `using CairoMakie`.
"""
function plot_cover end

"""
    plot_dhw(dom::SyntheticDomain, label, scenario)

DHW of every reef and the regional signal for one climate label and scenario, next to the
spread over the stochastic scenarios of every climate level. Requires `using CairoMakie`.
"""
function plot_dhw end

"""
    plot_cyclones(dom::SyntheticDomain, scenario)

Cyclone years and intensities for one scenario, and the mortality per functional group for
an average and the strongest event. Requires `using CairoMakie`.
"""
function plot_cyclones end

"""
    save_figures(dom::SyntheticDomain, dir)::Vector{String}

Save all plots as PNGs in `dir` and return their file names. Requires `using CairoMakie`.
"""
function save_figures end

"""
    plot_cscape_connectivity(dom::SyntheticDomain, options::CscapeOptions)

C~scape connectivity of the broadcast spawners and of the brooder, as `[destination,
source]` heatmaps. Requires `using CairoMakie`.
"""
function plot_cscape_connectivity end

"""
    plot_cscape_cyclones(dom::SyntheticDomain, options::CscapeOptions)

C~scape cyclone categories per year and replicate. Requires `using CairoMakie`.
"""
function plot_cscape_cyclones end

"""
    plot_cscape_cover(dom::SyntheticDomain, options::CscapeOptions)

Initial cover per reef for the six C~scape types (% of site area). Requires `using CairoMakie`.
"""
function plot_cscape_cover end

"""Whether the plotting extension is loaded (i.e. `using CairoMakie` has been run)."""
plots_available() = Base.get_extension(@__MODULE__, :SyntheticReefDomainCairoMakieExt) !== nothing

export DomainConfig, ReefConfig, ConnectivityConfig, CoverConfig, DHWConfig, CycloneConfig
export load_config, save_config, config_to_dict, validate
export SyntheticDomain, build_domain, write_domain, generate_domain
export build_reefs, build_connectivity, build_cover, build_dhw, build_waves, build_cyclones
export climate_labels, years, reef_outlines, FUNCTIONAL_GROUPS, ZERO_SCENARIO
export plot_map, plot_connectivity, plot_cover, plot_dhw, plot_cyclones, save_figures
export CscapeOptions, cscape_inputs, write_cscape_project, load_cscape_options,
    save_cscape_options, rscript_status, cyclone_category, CSCAPE_TYPES, CSCAPE_REGIONS
export domain_zip, cscape_zip
export plot_cscape_connectivity, plot_cscape_cyclones, plot_cscape_cover

end
