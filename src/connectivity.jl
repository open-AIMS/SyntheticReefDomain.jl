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

"""
    build_connectivity(cfg::DomainConfig, reefs::DataFrame)::Matrix{Float64}

Source x sink connectivity matrix: `self_recruitment` on the diagonal, `strength` times the
topology's link weights elsewhere, with each row's exports scaled down where needed so a
reef never exports more than `1 - self_recruitment` of its larvae.
"""
function build_connectivity(cfg::DomainConfig, reefs::DataFrame)::Matrix{Float64}
    c = cfg.connectivity
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
