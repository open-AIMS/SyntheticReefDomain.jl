# Connectivity: a source (row) x sink (column) larval transfer probability matrix.
#
# Note on how ADRIA uses it: ADRIA's own coral model normalises the matrix per sink column
# (`conn ./ sum(conn; dims=1)`, ADRIA/src/scenario.jl), i.e. only the *relative* makeup of
# each reef's incoming larval supply (home vs. each other reef) drives recruitment there.
# The absolute values still matter for ADRIA's connectivity-based decision criteria and
# for engines that move larvae through the matrix directly (CoralPop, Kora, C~scape).

const CONN_CUTOFF = 1e-6  # ADRIA zeroes values below this on load

"""Great-circle distance in km between two (lon, lat) points."""
function haversine_km(p::Tuple{Float64,Float64}, q::Tuple{Float64,Float64})::Float64
    R = 6371.0088
    dlat = deg2rad(q[2] - p[2])
    dlon = deg2rad(q[1] - p[1])
    a = sin(dlat / 2)^2 + cos(deg2rad(p[2])) * cos(deg2rad(q[2])) * sin(dlon / 2)^2
    return 2R * asin(min(1.0, sqrt(a)))
end

"""Pairwise distances (km) between reef centres."""
function distance_matrix_km(reefs::DataFrame)::Matrix{Float64}
    pts = collect(zip(reefs.lon, reefs.lat))
    return [haversine_km(p, q) for p in pts, q in pts]
end

# Relative link weights (0-1, zero diagonal) for each topology.
function _link_weights(c::ConnectivityConfig, reefs::DataFrame, rng)::Matrix{Float64}
    n = nrow(reefs)
    W = zeros(n, n)
    n == 1 && return W

    if c.topology in ("chain", "ring")
        wrap = c.topology == "ring" && n >= 3
        for i in 1:n, step in 1:c.reach
            w = 0.5^(step - 1)
            for (j, weight) in ((i + step, w), (i - step, w * (1 - c.asymmetry)))
                if 1 <= j <= n
                    W[i, j] = max(W[i, j], weight)
                elseif wrap
                    jj = mod1(j, n)
                    jj == i || (W[i, jj] = max(W[i, jj], weight))
                end
            end
        end
    elseif c.topology == "all_to_all"
        W .= 1.0
    elseif c.topology == "distance_decay"
        W .= exp.(-distance_matrix_km(reefs) ./ c.decay_km)
    elseif c.topology == "random"
        # Draw both numbers for every pair so link presence and weights stay stable when
        # `density` changes.
        for j in 1:n, i in 1:n
            present, weight = rand(rng), rand(rng)
            W[i, j] = present < c.density ? weight : 0.0
        end
    end  # "isolated": no links

    W[CartesianIndex.(1:n, 1:n)] .= 0.0
    return W
end

"""Number of links between different reefs (positive off-diagonal entries)."""
offdiag_links(C::AbstractMatrix)::Int =
    count(C[i, j] > 0 for i in axes(C, 1), j in axes(C, 2) if i != j)

"""Integral from `a` to `b` of the kernel `exp(-α x)` truncated to zero for `x < x0`."""
function _kernel_integral(a::Float64, b::Float64, α::Float64, x0::Float64)::Float64
    lo, hi = max(a, x0), max(b, x0)
    return hi > lo ? (exp(-α * lo) - exp(-α * hi)) / α : 0.0
end

# Distance-based connectivity for C~scape ("truncated dilution" model), following
# `wrangle_connectivity_per_group` in the AIMS script 03c_create_connectivity_files_v2.R.
# Differences: distances are great-circle distances between reef centres (the script uses
# polygon centroids in EPSG:3577), the truncation also applies to links between reefs, and
# the matrix is scaled to `matrix_total` instead of 1.
function _truncated_dilution(c::ConnectivityConfig, reefs::DataFrame)::Matrix{Float64}
    n = nrow(reefs)
    α = -log(0.05) / (1000 * c.dispersal_km)  # kernel falls to 5% at dispersal_km
    x0 = c.threshold_m
    A = reefs.area .* reefs.k ./ 100          # available substrate area (m²)
    r = sqrt.(A ./ π)                         # radius of a circle with that area (m)
    d = 1000 .* distance_matrix_km(reefs)     # centre distances (m)

    # Local retention: kernel from the centre to the edge of the reef. Links: kernel from
    # the near edge of the destination to its centre, weighted by the destination's area.
    M = [
        i == j ? _kernel_integral(0.0, r[i], α, x0) * A[i] :
        _kernel_integral(d[i, j] - r[j], d[i, j], α, x0) * A[j]
        for i in 1:n, j in 1:n
    ]

    total = sum(M)
    total > 0 && (M .*= c.matrix_total / total)
    M[M .< CONN_CUTOFF] .= 0.0
    return M
end

"""
    build_connectivity(cfg::DomainConfig, reefs::DataFrame)::Matrix{Float64}

Source x sink connectivity matrix: `self_recruitment` on the diagonal, `strength` times the
topology's link weights elsewhere, with each row's exports scaled down where needed so a
reef never exports more than `1 - self_recruitment` of its larvae.

The `truncated_dilution` topology instead follows the distance-based method used for
C~scape. For reefs with substrate area `A = area * k / 100`, radius `r = sqrt(A / π)` and
centre distance `d` (m), and a kernel `exp(-α x)` that is zero within `threshold_m` of the
source (`α = -ln(0.05) / dispersal_km`, so it falls to 5% at `dispersal_km`):

- diagonal (local retention) = kernel integrated from 0 to `r[i]`, times `A[i]`; zero when
  `r[i] <= threshold_m`;
- source `i` -> sink `j` = kernel integrated from `d - r[j]` to `d`, times `A[j]`;
- the whole matrix is then scaled to sum to `matrix_total`.
"""
function build_connectivity(cfg::DomainConfig, reefs::DataFrame)::Matrix{Float64}
    c = cfg.connectivity
    c.topology == "truncated_dilution" && return _truncated_dilution(c, reefs)
    n = nrow(reefs)
    C = c.strength .* _link_weights(c, reefs, component_rng(cfg.seed, "connectivity"))

    budget = 1.0 - c.self_recruitment
    for i in 1:n
        exports = sum(C[i, :])
        if exports > budget
            C[i, :] .*= exports > 0 ? budget / exports : 0.0
        end
        C[i, i] = c.self_recruitment
    end

    C[C .< CONN_CUTOFF] .= 0.0
    return clamp.(C, 0.0, 1.0)
end
