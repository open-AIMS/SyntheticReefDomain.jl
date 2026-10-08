# File writers. Each mirrors the structure of the matching file in
# ADRIA.jl/ADRIA/test/data/Test_domain, which is what ADRIA's loaders are tested against.
#
# NetCDF dimension order below is Julia's (the reverse of ncdump's), as in ADRIA.

const GENERATOR = "SyntheticReefDomain.jl"

function _gatts(title::AbstractString)
    return Dict{String,Any}(
        "title" => title,
        "source" => "$GENERATOR (synthetic data, not observations)",
        "history" => "created $(Dates.format(now(), "yyyy-mm-dd HH:MM"))"
    )
end

"""
    write_gpkg(path, reefs::DataFrame)

Reef polygons and attributes as a single-layer (`data`) EPSG:4326 geopackage. `k` is
written as a percentage (as in real ADRIA domains; ADRIA divides by 100 on load).
"""
function write_gpkg(path::AbstractString, reefs::DataFrame)
    GDF.write(
        path, reefs;
        layer_name="data", crs=GFT.EPSG(4326), geom_columns=(:geometry,)
    )
    return path
end

"""
    connectivity_filename(cfg)

ADRIA groups connectivity CSVs by the second `_`-separated token of the file name (a year)
and averages them, so the start year goes in that position.
"""
connectivity_filename(cfg::DomainConfig) = "connectivity_$(cfg.timeframe[1])_synthetic.csv"

"""
    write_connectivity(path, ids, C, cfg)

Source x sink matrix as CSV: `#` comment lines, a header whose first cell ADRIA ignores,
then one row per source reef. Values are written at full precision.
"""
function write_connectivity(
    path::AbstractString, ids::AbstractVector{<:AbstractString}, C::Matrix{Float64},
    cfg::DomainConfig
)
    c = cfg.connectivity
    open(path, "w") do io
        println(io, "# Synthetic transfer probability matrix ($GENERATOR).")
        println(io, "# Row names are the SOURCE reefs. Column names indicate the RECEIVING reef.")
        println(io, "# topology=$(c.topology) self_recruitment=$(c.self_recruitment) " *
                    "strength=$(c.strength) asymmetry=$(c.asymmetry) reach=$(c.reach) " *
                    "decay_km=$(c.decay_km) density=$(c.density)")
        println(io, join(["source_site"; ids], ","))
        for (i, id) in enumerate(ids)
            println(io, join([id; string.(C[i, :])], ","))
        end
    end
    return path
end

"""
    write_env_nc(path, varname, data, reefs; long_name, units)

A `timesteps x sites x member` environmental layer (DHW or wave Ub), laid out like the
Test_domain files: no coordinate variables, location IDs in a `reef_siteid` string
variable. ADRIA's loader finds the location dimension by matching its length to the IDs.
"""
function write_env_nc(
    path::AbstractString, varname::AbstractString, data::Array{Float64,3}, reefs::DataFrame;
    long_name::AbstractString, units::AbstractString
)
    T, n, S = size(data)
    d_time = NetCDF.NcDim("timesteps", T)
    d_site = NetCDF.NcDim("sites", n)
    d_member = NetCDF.NcDim("member", S)
    vars = [
        NetCDF.NcVar(
            varname, [d_time, d_site, d_member];
            t=Float64, atts=Dict{String,Any}("long_name" => long_name, "units" => units)
        ),
        NetCDF.NcVar("reef_siteid", [d_site]; t=String),
        NetCDF.NcVar("UNIQUE_ID", [d_site]; t=String),
        NetCDF.NcVar("longitude", [d_site]; t=Float64),
        NetCDF.NcVar("latitude", [d_site]; t=Float64)
    ]
    NetCDF.create(path, vars; gatts=_gatts(long_name)) do nc
        NetCDF.putvar(nc[varname], data)
        NetCDF.putvar(nc["reef_siteid"], collect(String, reefs.reef_siteid))
        NetCDF.putvar(nc["UNIQUE_ID"], collect(String, reefs.UNIQUE_ID))
        NetCDF.putvar(nc["longitude"], collect(Float64, reefs.lon))
        NetCDF.putvar(nc["latitude"], collect(Float64, reefs.lat))
    end
    return path
end

"""
    write_cover_nc(path, ids, cover)

Initial cover (`species x locations`, fraction of habitable area) with string coordinate
variables, which ADRIA needs to match functional groups by name.
"""
function write_cover_nc(path::AbstractString, ids::AbstractVector{<:AbstractString}, cover::Matrix{Float64})
    d_species = NetCDF.NcDim("species", N_GROUPS; values=copy(FUNCTIONAL_GROUPS))
    d_loc = NetCDF.NcDim("locations", length(ids); values=collect(String, ids))
    v = NetCDF.NcVar(
        "layer", [d_species, d_loc];
        t=Float64,
        atts=Dict{String,Any}(
            "long_name" => "Initial coral cover per functional group",
            "units" => "fraction of habitable area"
        )
    )
    NetCDF.create(path, [v]; gatts=_gatts("Initial coral cover (synthetic)")) do nc
        NetCDF.putvar(nc["layer"], cover)
    end
    return path
end

"""
    write_cyclone_nc(path, ids, mortality)

Cyclone mortality (`timesteps x locations x species x scenarios`) as a single data variable
with coordinate variables, as ADRIA reads it with `YAXArrays.Cube`.
"""
function write_cyclone_nc(
    path::AbstractString, ids::AbstractVector{<:AbstractString}, mortality::Array{Float64,4}
)
    T, n, G, S = size(mortality)
    dims = [
        NetCDF.NcDim("timesteps", T; values=collect(1:T)),
        NetCDF.NcDim("locations", n; values=collect(String, ids)),
        NetCDF.NcDim("species", G; values=copy(FUNCTIONAL_GROUPS)),
        NetCDF.NcDim("scenarios", S; values=collect(1:S))
    ]
    v = NetCDF.NcVar(
        "layer", dims;
        t=Float64,
        atts=Dict{String,Any}(
            "long_name" => "Cyclone mortality (system-wide events)",
            "units" => "fraction of cover lost"
        )
    )
    NetCDF.create(path, [v]; gatts=_gatts("Cyclone mortality (synthetic)")) do nc
        NetCDF.putvar(nc["layer"], mortality)
    end
    return path
end

"""
    write_datapackage(path, name, dom)

`datapackage.json` in the layout of ADRIA's Test_domain (version $(DPKG_VERSION)). ADRIA reads
`version`, `simulation_metadata.timeframe` and the column names of `spatial_data`; the
other entries document the package.
"""
function write_datapackage(path::AbstractString, name::AbstractString, dom::SyntheticDomain)
    cfg = dom.config
    labels = climate_labels(cfg)
    resources = Any[
        JSON.Object(
            "name" => "spatial_data",
            "description" => "Synthetic reef polygons and attributes",
            "path" => "spatial/$name.gpkg",
            "format" => "geopackage",
            "location_id_col" => "reef_siteid",
            "cluster_id_col" => "cluster_id",
            "k_col" => "k",
            "area_col" => "area"
        ),
        JSON.Object(
            "name" => "coral_cover",
            "description" => "Initial coral cover per functional group (fraction of habitable area)",
            "data" => ["locations", "species"],
            "path" => "spatial/coral_cover.nc",
            "format" => "netCDF"
        ),
        JSON.Object(
            "name" => "connectivity",
            "description" => "Synthetic connectivity ($(cfg.connectivity.topology)); rows are sources",
            "data" => ["row_source", "col_receiving"],
            "path" => "connectivity",
            "format" => "CSV"
        ),
        JSON.Object(
            "name" => "DHWs",
            "description" => "Synthetic degree heating weeks for climate label(s) $(join(labels, ", "))",
            "data" => ["reef_siteid", "dhw"],
            "path" => "DHWs",
            "format" => "netCDF"
        ),
        JSON.Object(
            "name" => "waves",
            "description" => "Wave Ub, constant per reef (no wave disturbance)",
            "data" => ["reef_siteid", "Ub"],
            "path" => "waves",
            "format" => "netCDF"
        ),
        JSON.Object(
            "name" => "cyclone_mortality",
            "description" => "System-wide synthetic cyclone mortality",
            "path" => "cyclones/cyclone_mortality.nc",
            "format" => "netCDF"
        )
    ]

    pkg = JSON.Object(
        "name" => name,
        "title" => "Synthetic reef system: $(cfg.reefs.n) reefs",
        "description" => "Synthetic ADRIA domain generated by $GENERATOR from " *
                         "generator_config.toml (seed $(cfg.seed)). Not observational data.",
        "version" => DPKG_VERSION,
        "sources" => Any[],
        "simulation_metadata" => JSON.Object("timeframe" => collect(cfg.timeframe)),
        "contributors" => Any[JSON.Object(
            "title" => GENERATOR,
            "role" => "author",
            "description" => "Generated synthetic data"
        )],
        "resources" => resources
    )
    open(path, "w") do io
        JSON.json(io, pkg; pretty=2)
    end
    return path
end

