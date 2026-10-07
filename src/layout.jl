# Reef positions, polygons and per-reef attributes.
#
# Positions are built in a local planar frame (x = km east, y = km north of the centre) and
# converted to lon/lat with an equirectangular approximation, which is accurate to well
# under 1% over the few tens of km a 20-reef system spans.

const KM_PER_DEG_LAT = 110.574
const KM_PER_DEG_LON_EQUATOR = 111.320

_to_lonlat(x_km, y_km, lon0, lat0) = (
    lon0 + x_km / (KM_PER_DEG_LON_EQUATOR * cosd(lat0)),
    lat0 + y_km / KM_PER_DEG_LAT
)
_to_km(lon, lat, lon0, lat0) = (
    (lon - lon0) * KM_PER_DEG_LON_EQUATOR * cosd(lat0),
    (lat - lat0) * KM_PER_DEG_LAT
)

# Unit vectors along the layout axis (bearing clockwise from north) and across it.
_axis(bearing_deg) = (sind(bearing_deg), cosd(bearing_deg))
_cross_axis(bearing_deg) = (-cosd(bearing_deg), sind(bearing_deg))

function _ring_positions(m::Int, spacing::Float64, bearing::Float64)
    m == 1 && return [(0.0, 0.0)]
    radius = spacing / (2 * sin(pi / m))
    return [
        (radius * sind(bearing + 360 * (i - 1) / m), radius * cosd(bearing + 360 * (i - 1) / m))
        for i in 1:m
    ]
end

"""Contiguous block sizes when splitting `n` reefs into `k` groups (sizes differ by <= 1)."""
_block_sizes(n::Int, k::Int) = [length(((i - 1) * n ÷ k + 1):(i * n ÷ k)) for i in 1:k]

"""
    reef_positions(r::ReefConfig, rng)::Tuple{Vector{Tuple{Float64,Float64}},Vector{Int}}

Reef centres in km relative to the system centre, in reef order, and the cluster index of
each reef. Reef order follows the layout so neighbouring indices are spatial neighbours.
"""
function reef_positions(r::ReefConfig, rng)
    n, s, b = r.n, r.spacing_km, r.bearing_deg
    u, v = _axis(b), _cross_axis(b)
    along(t) = (t * u[1], t * u[2])

    pos = if r.layout == "line"
        [along((i - (n + 1) / 2) * s) for i in 1:n]
    elseif r.layout == "grid"
        ncol = ceil(Int, sqrt(n))
        nrow = ceil(Int, n / ncol)
        cells = Tuple{Float64,Float64}[]
        for row in 1:nrow, j in 1:ncol
            # Serpentine order keeps consecutive reefs adjacent (useful for chain topologies)
            col = isodd(row) ? j : ncol - j + 1
            a = (col - (ncol + 1) / 2) * s
            c = (row - (nrow + 1) / 2) * s
            push!(cells, (a * u[1] + c * v[1], a * u[2] + c * v[2]))
        end
        cells[1:n]
    elseif r.layout == "ring"
        _ring_positions(n, s, b)
    elseif r.layout == "random"
        side = s * sqrt(n) * 1.5
        min_dist = 0.6 * s
        pts = Tuple{Float64,Float64}[]
        for _ in 1:n
            local p
            for _ in 1:1000
                p = ((rand(rng) - 0.5) * side, (rand(rng) - 0.5) * side)
                all(q -> hypot(p[1] - q[1], p[2] - q[2]) >= min_dist, pts) && break
            end
            push!(pts, p)
        end
        # Order along the main axis so reef order (and gradients) runs across the system
        sort!(pts; by=p -> p[1] * u[1] + p[2] * u[2])
    else  # clustered
        k = r.n_clusters
        sizes = _block_sizes(n, k)
        pts = Tuple{Float64,Float64}[]
        for (j, m) in enumerate(sizes)
            cx, cy = along((j - (k + 1) / 2) * r.cluster_spacing_km)
            for (px, py) in _ring_positions(m, s, b + 360 * rand(rng))
                jitter = 0.15 * s
                push!(pts, (cx + px + jitter * randn(rng), cy + py + jitter * randn(rng)))
            end
        end
        pts
    end

    sizes = _block_sizes(n, r.n_clusters)
    cluster = reduce(vcat, [fill(j, m) for (j, m) in enumerate(sizes)])
    return pos, cluster
end

"""
    reef_polygon(center_km, area_m2, rng; n_vertices=40)

An irregular, roughly elliptical outline with exactly `area_m2`, as km offsets.
"""
function reef_polygon(center_km::Tuple{Float64,Float64}, area_m2::Float64, rng; n_vertices=40)
    aspect = 1.2 + 0.8 * rand(rng)
    rot = 2pi * rand(rng)
    ph1, ph2 = 2pi * rand(rng), 2pi * rand(rng)
    angles = range(0, 2pi; length=n_vertices + 1)[1:n_vertices]

    pts = map(angles) do t
        radius = 1 + 0.10 * sin(2t + ph1) + 0.06 * sin(3t + ph2)
        x, y = radius * aspect * cos(t), radius * sin(t) / aspect
        (x * cos(rot) - y * sin(rot), x * sin(rot) + y * cos(rot))
    end

    # Shoelace area, then scale to the requested area (km^2)
    a = 0.5 * abs(sum(
        pts[i][1] * pts[mod1(i + 1, end)][2] - pts[mod1(i + 1, end)][1] * pts[i][2]
        for i in eachindex(pts)
    ))
    scale = sqrt((area_m2 / 1e6) / a)
    return [(center_km[1] + scale * p[1], center_km[2] + scale * p[2]) for p in pts]
end

function _distinct_pairs!(ub::Vector{Float64}, depth::Vector{Float64})
    seen = Set{Tuple{Float64,Float64}}()
    for i in eachindex(ub)
        while (ub[i], depth[i]) in seen
            ub[i] += 1e-4
        end
        push!(seen, (ub[i], depth[i]))
    end
    return ub
end

"""
    reef_outlines(reefs::DataFrame, center; units=:lonlat)

Polygon vertices of each reef as `(lon, lat)` tuples, or `(x, y)` km offsets from `center`
with `units=:km` (convenient for plotting with equal axis scaling).
"""
function reef_outlines(reefs::DataFrame, center; units::Symbol=:lonlat)
    units in (:lonlat, :km) || throw(ArgumentError("units must be :lonlat or :km"))
    lon0, lat0 = center
    return map(reefs.geometry) do geom
        ring = AG.getgeom(geom, 0)
        pts = [(AG.getx(ring, i), AG.gety(ring, i)) for i in 0:(AG.ngeom(ring) - 1)]
        units == :km ? [_to_km(lon, lat, lon0, lat0) for (lon, lat) in pts] : pts
    end
end

"""
    reef_id(i::Int)::String

Location identifier of reef `i`. Zero padded (`reef_01`) because ADRIA sorts locations by
their string ID.
"""
reef_id(i::Int)::String = @sprintf("reef_%02d", i)

"""
    build_reefs(cfg::DomainConfig)::DataFrame

One row per reef with its position, polygon (`geometry`, EPSG:4326) and the attributes
ADRIA and C~scape read (`reef_siteid`, `UNIQUE_ID`, `cluster_id`, `k`, `area`, `depth_med`,
`ub_med`, `CB_CALIB_GROUPS`) plus informational columns (`habitat`, DHW exposure, ...).
"""
function build_reefs(cfg::DomainConfig)::DataFrame
    r = cfg.reefs
    n = r.n
    rng(tag) = component_rng(cfg.seed, tag)

    pos, cluster = reef_positions(r, rng("layout"))
    area = resolve_attr(r.area_m2, n, rng("area_m2"))
    k = resolve_attr(r.k_percent, n, rng("k_percent"))
    depth = resolve_attr(r.depth_med, n, rng("depth_med"))
    ub = resolve_attr(r.ub_med, n, rng("ub_med"))
    r.distinct_ub_depth && _distinct_pairs!(ub, depth)
    exposure = resolve_attr(cfg.dhw.reef_exposure, n, rng("dhw_exposure"))
    offset = resolve_attr(cfg.dhw.reef_offset, n, rng("dhw_offset"))

    lon0, lat0 = r.center
    poly_rng = rng("polygons")
    outlines = [reef_polygon(pos[i], area[i], poly_rng) for i in 1:n]
    geometry = map(outlines) do outline
        ring = [_to_lonlat(x, y, lon0, lat0) for (x, y) in outline]
        push!(ring, first(ring))  # close the ring
        AG.createpolygon(ring)
    end
    lonlat = [_to_lonlat(p[1], p[2], lon0, lat0) for p in pos]
    ids = reef_id.(1:n)

    return DataFrame(
        reef_siteid=ids,
        UNIQUE_ID=ids,
        site_id=collect(1:n),
        cluster_id=[@sprintf("cluster_%02d", c) for c in cluster],
        habitat=resolve_categorical(r.habitat, n, String),
        zone_type=fill("General Use", n),
        k=k,
        area=area,
        depth_med=depth,
        depth_mean=copy(depth),
        ub_med=ub,
        CB_CALIB_GROUPS=resolve_categorical(r.cb_calib_group, n, Int),
        dhw_exposure=exposure,
        dhw_offset=offset,
        lon=first.(lonlat),
        lat=last.(lonlat),
        x_km=first.(pos),
        y_km=last.(pos),
        geometry=geometry
    )
end
