# Initial coral cover per functional group, as a fraction of each reef's habitable area
# (ADRIA's native unit; ADRIA splits each group into its 7 size classes itself).

"""
    build_cover(cfg::DomainConfig, reefs::DataFrame)::Matrix{Float64}

Initial cover as a `5 x n` matrix (groups in `FUNCTIONAL_GROUPS` order x reefs), in
fractions of habitable area, so each column sums to at most 1.
"""
function build_cover(cfg::DomainConfig, reefs::DataFrame)::Matrix{Float64}
    cc = cfg.cover
    n = nrow(reefs)

    total = resolve_attr(cc.total_percent, n, component_rng(cfg.seed, "cover_total")) ./ 100
    if cc.basis == "reef"
        k = reefs.k ./ 100
        too_high = findall(total .> k .+ 1e-12)
        isempty(too_high) || throw(ArgumentError(
            "cover.total_percent exceeds the habitable fraction (k) on " *
            join(["$(reefs.reef_siteid[i]) ($(round(100total[i]; digits=1))% > k = " *
                  "$(round(100k[i]; digits=1))%)" for i in too_high], ", ") *
            ". With basis = \"reef\" cover is a percentage of the whole reef area."
        ))
        total = total ./ k
    end

    w = cc.composition ./ sum(cc.composition)
    composition = repeat(w, 1, n)
    if cc.jitter > 0
        # Vary the mix per reef around `w`; groups with zero weight stay at zero.
        rng = component_rng(cfg.seed, "cover_composition")
        active = w .> 0
        if count(active) > 1
            dist = Dirichlet(w[active] ./ cc.jitter)
            for j in 1:n
                composition[active, j] .= rand(rng, dist)
            end
        end
    end

    return composition .* reshape(total, 1, n)
end
