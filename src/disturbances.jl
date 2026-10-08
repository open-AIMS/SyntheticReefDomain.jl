# Disturbance time series: heat stress (DHW), cyclone mortality and (static) wave exposure.
# Arrays use ADRIA's (Julia) dimension order: timesteps x locations [x groups] x scenarios.
#
# Scenario 1 of every DHW and cyclone array is disturbance-free (all zeros): a control
# scenario present in every domain. The configured stochastic scenarios follow as 2..S+1,
# so the files hold `n_dhw_scenarios + 1` and `n_cyclone_scenarios + 1` scenarios.

"""Index of the disturbance-free scenario in every DHW, wave and cyclone file."""
const ZERO_SCENARIO = 1

"""Add the all-zero scenario in front of the last (scenario) dimension."""
function _prepend_zero_scenario(A::AbstractArray{Float64})
    return cat(zeros(size(A)[1:(end - 1)]..., 1), A; dims=ndims(A))
end

"""
    build_dhw(cfg::DomainConfig, reefs::DataFrame)

DHW for every climate level. Returns `(dhw, regional)`: dictionaries keyed by climate label
holding the `T x n x (S + 1)` reef DHW and the `T x (S + 1)` regional signal it was built
from, where `S = n_dhw_scenarios`. Scenario 1 is all zeros (no heat stress); the stochastic
scenarios are 2..S+1.

All climate levels use the same underlying random draws (common random numbers), and each
reef's local noise comes from its own stream, so adding reefs or climate levels leaves
existing series unchanged.
"""
function build_dhw(cfg::DomainConfig, reefs::DataFrame)
    yrs = years(cfg)
    T, n, S = length(yrs), nrow(reefs), cfg.n_dhw_scenarios

    rng = component_rng(cfg.seed, "dhw_regional")
    z = randn(rng, T, S)          # AR(1) innovations
    u_event = rand(rng, T, S)     # event occurrence
    u_mag = rand(rng, T, S)       # event magnitude
    local_noise = Array{Float64}(undef, T, n, S)
    for r in 1:n
        local_noise[:, r, :] .= randn(component_rng(cfg.seed, "dhw_local_$(r)"), T, S)
    end

    dhw = Dict{Int,Array{Float64,3}}()
    regional = Dict{Int,Matrix{Float64}}()
    for label in climate_labels(cfg)
        p = climate_params(cfg.dhw, label)
        innovation_sd = p.interannual_sd * sqrt(1 - p.ar1^2)
        prescribed = in.(yrs, Ref(Set(p.event_years)))
        lo, hi = p.event_magnitude

        R = Matrix{Float64}(undef, T, S)
        for s in 1:S
            anomaly = p.interannual_sd * z[1, s]
            for t in 1:T
                t > 1 && (anomaly = p.ar1 * anomaly + innovation_sd * z[t, s])
                is_event = prescribed[t] || u_event[t, s] < p.event_prob
                event = is_event ? lo + (hi - lo) * u_mag[t, s] : 0.0
                R[t, s] = p.baseline + p.trend_per_decade * (t - 1) / 10 + anomaly + event
            end
        end

        D = Array{Float64}(undef, T, n, S)
        for s in 1:S, r in 1:n, t in 1:T
            D[t, r, s] = max(
                0.0,
                reefs.dhw_exposure[r] * R[t, s] + reefs.dhw_offset[r] +
                p.local_sd * local_noise[t, r, s]
            )
        end

        dhw[label] = _prepend_zero_scenario(D)
        regional[label] = _prepend_zero_scenario(R)
    end

    return dhw, regional
end

"""
    build_cyclones(cfg::DomainConfig, n_reefs::Int)

System-wide cyclone mortality. Returns `(mortality, intensity)`: the `T x n x 5 x (S + 1)`
mortality fraction ADRIA reads, identical on every reef, and the `T x (S + 1)` event
intensity (0 in years without a cyclone), where `S = n_cyclone_scenarios`. Scenario 1 has no
cyclones; `event_years` apply to the stochastic scenarios 2..S+1.
"""
function build_cyclones(cfg::DomainConfig, n_reefs::Int)
    cy = cfg.cyclones
    yrs = years(cfg)
    T, S = length(yrs), cfg.n_cyclone_scenarios

    rng = component_rng(cfg.seed, "cyclones")
    u_event = rand(rng, T, S)
    u_intensity = rand(rng, T, S)

    intensity = zeros(T, S)
    if cy.enabled
        prescribed = in.(yrs, Ref(Set(cy.event_years)))
        lo, hi = cy.intensity
        for s in 1:S, t in 1:T
            if prescribed[t] || u_event[t, s] < cy.annual_prob
                intensity[t, s] = lo + (hi - lo) * u_intensity[t, s]
            end
        end
    end

    mortality = zeros(T, n_reefs, N_GROUPS, S)
    for s in 1:S, g in 1:N_GROUPS, t in 1:T
        intensity[t, s] > 0 || continue
        mortality[t, :, g, s] .= clamp(intensity[t, s] * cy.group_vulnerability[g], 0.0, 1.0)
    end

    return _prepend_zero_scenario(mortality), _prepend_zero_scenario(intensity)
end

"""
    build_waves(cfg::DomainConfig, reefs::DataFrame)::Array{Float64,3}

Wave Ub (`T x n x (S + 1)`, with `S = n_dhw_scenarios`). There is no wave disturbance: each
reef's `ub_med` is repeated over time and scenarios (or zeros with `waves_from_ub_med = false`),
including the disturbance-free scenario 1.

A file is written anyway because ADRIA indexes waves with the DHW scenario index during
site selection, so the wave file must have as many scenarios as the DHW file.
"""
function build_waves(cfg::DomainConfig, reefs::DataFrame)::Array{Float64,3}
    T, n, S = length(years(cfg)), nrow(reefs), cfg.n_dhw_scenarios + 1
    ub = cfg.reefs.waves_from_ub_med ? reefs.ub_med : zeros(n)
    return repeat(reshape(ub, 1, n, 1), T, 1, S)
end
