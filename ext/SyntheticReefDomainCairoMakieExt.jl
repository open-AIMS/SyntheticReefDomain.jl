# Plots of a synthetic domain (the same figures the Pluto notebook shows), and saving them
# as PNGs for the domain README. Loaded automatically with `using CairoMakie`.
module SyntheticReefDomainCairoMakieExt

using CairoMakie
using Statistics: mean, quantile
using DataFrames: nrow, names

import SyntheticReefDomain as SRD
using SyntheticReefDomain: SyntheticDomain, ZERO_SCENARIO, climate_labels, years, reef_outlines

const GROUP_LABELS = [
    "Tabular Acropora", "Corymbose Acropora", "Corymbose non-Acropora",
    "Small massives", "Large massives"
]

short_ids(dom::SyntheticDomain) = last.(split.(dom.reefs.reef_siteid, "_"))

function SRD.plot_map(dom::SyntheticDomain)
    reefs = dom.reefs
    C = dom.connectivity
    n = size(C, 1)
    x, y = reefs.x_km, reefs.y_km
    cover = vec(sum(dom.cover; dims=1)) .* 100

    # Map extent with a usable aspect ratio for elongated layouts (e.g. a line of reefs);
    # the axis is sized in pixels so 1 km is equally long on both axes.
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    pad = 0.12 * max(xmax - xmin, ymax - ymin, 2 * dom.config.reefs.spacing_km)
    w = max(xmax - xmin, 0.3 * (ymax - ymin)) + 2pad
    h = max(ymax - ymin, 0.3 * (xmax - xmin)) + 2pad
    cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
    px = 620
    ax_w, ax_h = w >= h ? (px, px * h / w) : (px * w / h, px)

    fig = Figure()
    ax = Axis(fig[1, 1];
        width=ax_w, height=ax_h,
        xlabel="km east of centre", ylabel="km north of centre",
        title="Reefs and connectivity",
        subtitle="colour = initial cover, size = area, outline = cluster"
    )
    limits!(ax, cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2)

    # Connectivity as curved arrows (so i -> j and j -> i do not overlap)
    off = [C[i, j] for i in 1:n, j in 1:n if i != j]
    cmax = isempty(off) ? 0.0 : maximum(off)
    if cmax > 0
        for i in 1:n, j in 1:n
            (i == j || C[i, j] < 0.05 * cmax) && continue
            w_ij = C[i, j] / cmax
            p0, p1 = Point2(x[i], y[i]), Point2(x[j], y[j])
            d = p1 - p0
            mid = (p0 + p1) / 2 + 0.15 * Point2(-d[2], d[1])
            bez(t) = (1 - t)^2 * p0 + 2 * (1 - t) * t * mid + t^2 * p1
            lines!(ax, [bez(t) for t in range(0, 1; length=24)];
                color=(:steelblue, 0.2 + 0.6w_ij), linewidth=0.5 + 3w_ij)
            tangent = 2 * 0.25 * (mid - p0) + 2 * 0.75 * (p1 - mid)
            scatter!(ax, [bez(0.75)];
                marker=:utriangle, markersize=7 + 8w_ij, color=(:steelblue, 0.35 + 0.65w_ij),
                rotation=atan(tangent[2], tangent[1]) - pi / 2)
        end
    end

    # Reef polygons (true size; often tiny next to the spacing) and markers
    for outline in reef_outlines(dom; units=:km)
        poly!(ax, Point2.(outline); color=(:black, 0.25), strokewidth=0)
    end
    clusters = sort(unique(reefs.cluster_id))
    stroke = [Makie.wong_colors()[mod1(findfirst(==(c), clusters), 7)] for c in reefs.cluster_id]
    size_px = 12 .+ 14 .* sqrt.(reefs.area ./ maximum(reefs.area))
    sc = scatter!(ax, x, y;
        color=cover, colormap=:viridis, colorrange=(0, max(1, maximum(cover))),
        markersize=size_px, strokecolor=stroke, strokewidth=2.5)
    text!(ax, x, y; text=short_ids(dom), align=(:center, :bottom), offset=(0, 12), fontsize=11)
    Colorbar(fig[1, 2], sc; label="initial cover (% of habitable area)", height=max(ax_h, 220))
    resize_to_layout!(fig)
    return fig
end

function SRD.plot_connectivity(dom::SyntheticDomain)
    C = dom.connectivity
    n = size(C, 1)
    ids = short_ids(dom)
    fig = Figure(; size=(460, 400))
    ax = Axis(fig[1, 1];
        xlabel="sink (receiving) reef", ylabel="source reef", title="Connectivity",
        xticks=(1:n, ids), yticks=(1:n, ids), yreversed=true)
    hm = heatmap!(ax, 1:n, 1:n, permutedims(C);
        colormap=:Blues, colorrange=(0, max(maximum(C), 1e-9)))
    Colorbar(fig[1, 2], hm; label="transfer probability")
    return fig
end

function SRD.plot_cover(dom::SyntheticDomain)
    n = size(dom.cover, 2)
    xs = repeat(1:n; inner=5)
    groups = repeat(1:5; outer=n)
    colors = Makie.wong_colors()[1:5]
    fig = Figure(; size=(820, 320))
    ax = Axis(fig[1, 1];
        xticks=(1:n, short_ids(dom)), xlabel="reef",
        ylabel="% of habitable area", title="Initial cover by functional group")
    barplot!(ax, xs, vec(dom.cover) .* 100; stack=groups, color=colors[groups])
    Legend(fig[1, 2], [PolyElement(; color=c) for c in colors], GROUP_LABELS; framevisible=false)
    return fig
end

function SRD.plot_dhw(dom::SyntheticDomain, label::Integer, scen::Integer)
    yrs = collect(years(dom))
    D, R = dom.dhw[label], dom.dhw_regional[label]
    exposure = dom.reefs.dhw_exposure
    lo, hi = extrema(exposure)
    fig = Figure(; size=(940, 340))
    title = scen == ZERO_SCENARIO ? "Label $label, scenario 1: no heat stress (control)" :
        "Label $label, scenario $scen (colour = reef exposure)"
    ax1 = Axis(fig[1, 1]; title=title, xlabel="year", ylabel="DHW")
    for r in axes(D, 2)
        lines!(ax1, yrs, D[:, r, scen];
            color=exposure[r], colormap=:thermal, colorrange=(lo, hi + 1e-9), linewidth=1)
    end
    lines!(ax1, yrs, R[:, scen]; color=:black, linewidth=2.5, label="regional signal")
    axislegend(ax1; position=:lt)

    ax2 = Axis(fig[1, 2];
        title="Stochastic scenarios: mean over reefs (median, 10-90%)", xlabel="year")
    for (k, l) in enumerate(climate_labels(dom))
        stochastic = dom.dhw[l][:, :, (ZERO_SCENARIO + 1):end]
        m = dropdims(mean(stochastic; dims=2); dims=2)
        q = [quantile(m[t, :], p) for t in axes(m, 1), p in (0.1, 0.5, 0.9)]
        c = Makie.wong_colors()[mod1(k, 7)]
        band!(ax2, yrs, q[:, 1], q[:, 3]; color=(c, 0.25))
        lines!(ax2, yrs, q[:, 2]; color=c, linewidth=2, label="label $l")
    end
    axislegend(ax2; position=:lt)
    linkyaxes!(ax1, ax2)
    return fig
end

function SRD.plot_cyclones(dom::SyntheticDomain, scen::Integer)
    yrs = collect(years(dom))
    intensity = dom.cyclone_intensity
    cy = dom.config.cyclones
    fig = Figure(; size=(940, 320))
    n_events = vec(count(>(0), intensity; dims=1))
    mean_events = round(mean(n_events[(ZERO_SCENARIO + 1):end]); digits=1)
    title = scen == ZERO_SCENARIO ? "Scenario 1: no cyclones (control)" :
        "Scenario $scen: $(n_events[scen]) cyclone(s); mean $mean_events per stochastic scenario"
    ax1 = Axis(fig[1, 1];
        title=title, xlabel="year", ylabel="intensity", limits=(nothing, (0, 1.05)))
    barplot!(ax1, yrs, intensity[:, scen]; color=:indianred, gap=0)
    ax2 = Axis(fig[1, 2];
        title="Mortality per event", ylabel="fraction of cover lost",
        xticks=(1:5, ["Tab. Acro.", "Cor. Acro.", "Cor. non-Acro.", "Small mass.", "Large mass."]),
        xticklabelrotation=pi / 6, limits=(nothing, (0, 1.05)))
    lo, hi = cy.intensity
    barplot!(ax2, 1:5, clamp.(cy.group_vulnerability .* hi, 0, 1);
        color=(:indianred, 0.35), label="strongest event")
    barplot!(ax2, 1:5, clamp.(cy.group_vulnerability .* (lo + hi) / 2, 0, 1);
        color=:indianred, label="average event")
    axislegend(ax2; position=:rt)
    colsize!(fig.layout, 2, Relative(0.35))
    return fig
end

# --- C~scape export previews ------------------------------------------------------------

function SRD.plot_cscape_connectivity(dom::SyntheticDomain, options::SRD.CscapeOptions)
    C = SRD.cscape_inputs(dom, options).connectivity
    n = size(C, 1)
    ids = short_ids(dom)
    fig = Figure(; size=(900, 400))
    for (k, (t, title)) in enumerate(((1, "Broadcast spawners (5 types)"), (size(C, 3), "Brooder")))
        M = C[:, :, t]
        ax = Axis(fig[1, 2k - 1];
            title=title, xlabel="source reef", ylabel="destination reef",
            xticks=(1:n, ids), yticks=(1:n, ids), yreversed=true)
        hm = heatmap!(ax, 1:n, 1:n, permutedims(M);
            colormap=:Blues, colorrange=(0, max(maximum(M), 1e-12)))
        Colorbar(fig[1, 2k], hm; label="transfer probability")
    end
    return fig
end

function SRD.plot_cscape_cyclones(dom::SyntheticDomain, options::SRD.CscapeOptions)
    inputs = SRD.cscape_inputs(dom, options)
    cats = inputs.disturbances[first(climate_labels(dom))][:, 1, :, 3]  # years x replicates
    yrs = collect(years(dom))
    fig = Figure(; size=(940, 320))
    ax = Axis(fig[1, 1];
        title="Cyclone category per year and replicate (same on every reef; replicate 1 = control)",
        xlabel="year", ylabel="replicate (= DHW scenario)")
    hm = heatmap!(ax, yrs, 1:size(cats, 2), cats;
        colormap=cgrad([:white, :mistyrose, :salmon, :indianred, :firebrick, :darkred], 6; categorical=true),
        colorrange=(-0.5, 5.5))
    Colorbar(fig[1, 2], hm; label="category", ticks=0:5)
    return fig
end

function SRD.plot_cscape_cover(dom::SyntheticDomain, options::SRD.CscapeOptions)
    ic = SRD.cscape_inputs(dom, options).init_cover
    types = names(ic)[2:end]
    n, m = nrow(ic), length(types)
    values = Matrix(ic[:, 2:end])  # reefs x types
    xs = repeat(1:n; inner=m)
    groups = repeat(1:m; outer=n)
    colors = [Makie.wong_colors()[1:5]; [:gray60]]
    fig = Figure(; size=(820, 320))
    ax = Axis(fig[1, 1];
        xticks=(1:n, short_ids(dom)), xlabel="reef", ylabel="% of site area",
        title="C~scape initial cover by type (initialisationData.csv)")
    barplot!(ax, xs, vec(permutedims(values)); stack=groups, color=colors[groups])
    Legend(fig[1, 2], [PolyElement(; color=c) for c in colors], types; framevisible=false)
    return fig
end

function SRD.save_figures(dom::SyntheticDomain, dir::AbstractString)
    mkpath(dir)
    scen = ZERO_SCENARIO + 1  # first stochastic scenario
    figures = Tuple{String,Figure}[
        ("map.png", SRD.plot_map(dom)),
        ("connectivity.png", SRD.plot_connectivity(dom)),
        ("cover.png", SRD.plot_cover(dom))
    ]
    for label in climate_labels(dom)
        push!(figures, ("dhw_$(label).png", SRD.plot_dhw(dom, label, scen)))
    end
    push!(figures, ("cyclones.png", SRD.plot_cyclones(dom, scen)))

    for (file, fig) in figures
        save(joinpath(dir, file), fig; px_per_unit=2)
    end
    return first.(figures)
end

end
