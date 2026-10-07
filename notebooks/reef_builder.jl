### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 65ed2d60-c20a-11f1-203a-a7b44f43174d
md"""
# Synthetic reef builder

Build a synthetic ADRIA domain with 1-20 reefs. Pick a starting configuration, adjust the
reef system, connectivity, initial cover and disturbances, check the previews, then
**write the domain** -- or save the settings as a TOML config and regenerate exactly the
same domain from a script with `generate_domain("configs/<name>.toml")`.

Every change rebuilds the domain in memory; nothing is written to disk until you press one
of the buttons at the bottom.
"""

# ╔═╡ 65eedb10-c20a-11f1-049d-5511fec48719
begin
	import Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	using SyntheticReefDomain, CairoMakie, Statistics
	# Explicit imports: Makie also exports a `Slider`
	using PlutoUI: PlutoUI, Slider, Select, CheckBox, TextField, NumberField, CounterButton
	import SyntheticReefDomain as SRD
	CairoMakie.activate!(; type="png", px_per_unit=2)
	const REPO = normpath(joinpath(@__DIR__, ".."))
	const CONFIG_DIR = joinpath(REPO, "configs")
end

# ╔═╡ 65eedb10-c20a-11f1-0155-a1518672437e
config_files = let
	files = String[]
	for dir in (CONFIG_DIR, joinpath(CONFIG_DIR, "examples"))
		isdir(dir) || continue
		append!(files, [joinpath(dir, f) for f in readdir(dir) if endswith(f, ".toml")])
	end
	sort!([relpath(f, CONFIG_DIR) for f in files]; by=f -> (f != "default.toml", f))
end

# ╔═╡ 65eedb10-c20a-11f1-0fff-792e6ec4859a
md"""
## Starting point
Start from $(@bind base_file Select(config_files)) and
$(@bind keep_lists CheckBox(default=true)) keep per-reef lists from that file (used when
their length matches the number of reefs; the matching range controls below are then ignored).
"""

# ╔═╡ 65eedb10-c20a-11f1-37f7-652c033e5226
base = load_config(joinpath(CONFIG_DIR, base_file))

# ╔═╡ 65eedb10-c20a-11f1-3a88-459f1b87407e
begin
	"Closest value of `r` to `v` (widget defaults must lie on the widget's grid)."
	snap(v, r) = r[argmin(abs.(r .- v))]

	"A min-max pair of number fields initialised from a per-reef attribute specification."
	function range_input(spec, r)
		lo, hi = SRD.attr_bounds(spec)
		return PlutoUI.combine() do Child
			md"$(Child(:min, NumberField(r; default=snap(lo, r)))) to $(Child(:max, NumberField(r; default=snap(hi, r))))"
		end
	end

	"Per-reef specification from a min-max pair (a single value when min == max)."
	function range_spec(x)
		lo, hi = minmax(Float64(x.min), Float64(x.max))
		return lo == hi ? lo : (min=lo, max=hi)
	end

	parse_years(s) = [parse(Int, t) for t in split(strip(s), r"[,;\s]+") if !isempty(t)]

	"Parse `label = trend; label = trend`, keeping other overrides of existing levels."
	function parse_climate(s, base_climate)
		out = Dict{Int,Dict{String,Any}}()
		for item in split(s, r"[;\n]")
			isempty(strip(item)) && continue
			parts = strip.(split(item, "="))
			length(parts) == 2 || error("Climate levels must look like `45 = 1.2; 85 = 2.0`, got `$item`.")
			label = parse(Int, parts[1])
			overrides = copy(get(base_climate, label, Dict{String,Any}()))
			overrides["trend_per_decade"] = parse(Float64, parts[2])
			out[label] = overrides
		end
		return out
	end

	const HABITATS = ["Crest", "Outer Flat", "Sheltered Slope", "Slope"]
	const GROUP_LABELS = [
		"Tabular Acropora", "Corymbose Acropora", "Corymbose non-Acropora",
		"Small massives", "Large massives"
	]
	short_ids(dom) = last.(split.(dom.reefs.reef_siteid, "_"))
end

# ╔═╡ 65eedb10-c20a-11f1-37ed-37e36fc02c3c
@bind dom_ui PlutoUI.combine() do Child
	md"""
	## Domain
	- **Name** (also the folder name) $(Child(:name, TextField(default=base.name)))
	- **Years** $(Child(:y0, NumberField(1950:2300; default=base.timeframe[1]))) to $(Child(:y1, NumberField(1950:2300; default=base.timeframe[2]))) (at least 26 years)
	- **Seed** $(Child(:seed, NumberField(0:1_000_000_000; default=base.seed)))
	- **DHW scenarios** $(Child(:n_dhw, Slider(1:100; default=clamp(base.n_dhw_scenarios, 1, 100), show_value=true)))
	- **Cyclone scenarios** $(Child(:n_cyc, Slider(1:100; default=clamp(base.n_cyclone_scenarios, 1, 100), show_value=true)))
	"""
end

# ╔═╡ 65eedb10-c20a-11f1-3921-d9d49b19e4ae
@bind reef_ui PlutoUI.combine() do Child
	r = base.reefs
	habitat = r.habitat isa AbstractString && r.habitat in HABITATS ? r.habitat : "Slope"
	cb = r.cb_calib_group isa Integer ? r.cb_calib_group : 1
	md"""
	## Reefs
	- **Number of reefs** $(Child(:n, Slider(1:20; default=r.n, show_value=true)))
	- **Layout** $(Child(:layout, Select(collect(SRD.LAYOUTS); default=r.layout))) with $(Child(:n_clusters, Slider(1:20; default=r.n_clusters, show_value=true))) cluster(s)
	- **Spacing** $(Child(:spacing, Slider(0.5:0.5:50; default=snap(r.spacing_km, 0.5:0.5:50), show_value=true))) km; between cluster centres $(Child(:cluster_spacing, Slider(1:100; default=snap(r.cluster_spacing_km, 1:100), show_value=true))) km
	- **Bearing** of the main axis $(Child(:bearing, Slider(0:5:355; default=snap(mod(r.bearing_deg, 360), 0:5:355), show_value=true)))° from north
	- **Centre** lon $(Child(:lon, NumberField(-180:0.01:180; default=snap(r.center[1], -180:0.01:180)))) lat $(Child(:lat, NumberField(-80:0.01:80; default=snap(r.center[2], -80:0.01:80))))
	- **Area** (m²) $(Child(:area, range_input(r.area_m2, 1000:1000:5_000_000)))
	- **Habitable fraction k** (%) $(Child(:k, range_input(r.k_percent, 1:100)))
	- **Median depth** (m) $(Child(:depth, range_input(r.depth_med, 0.5:0.5:40)))
	- **Median wave orbital velocity `ub_med`** (m/s, read by C~scape) $(Child(:ub, range_input(r.ub_med, 0:0.01:2)))
	- **Habitat** $(Child(:habitat, Select(HABITATS; default=habitat))) **CoralBlox calibration group** $(Child(:cb, NumberField(1:12; default=cb)))
	- $(Child(:waves_from_ub, CheckBox(default=r.waves_from_ub_med))) write each reef's `ub_med` as its (constant) wave series
	"""
end

# ╔═╡ 65eedb10-c20a-11f1-3f17-1d1b42bb3704
@bind conn_ui PlutoUI.combine() do Child
	c = base.connectivity
	md"""
	## Connectivity
	Source (row) → sink (column) transfer probabilities.
	- **Topology** $(Child(:topology, Select(collect(SRD.TOPOLOGIES); default=c.topology)))
	- **Self-recruitment** $(Child(:self, Slider(0:0.01:1; default=snap(c.self_recruitment, 0:0.01:1), show_value=true))) **Strength** (largest off-diagonal value) $(Child(:strength, Slider(0:0.005:1; default=snap(c.strength, 0:0.005:1), show_value=true)))
	- *chain / ring:* **asymmetry** (0 symmetric, 1 downstream only) $(Child(:asymmetry, Slider(0:0.05:1; default=snap(c.asymmetry, 0:0.05:1), show_value=true))) **reach** $(Child(:reach, Slider(1:5; default=clamp(c.reach, 1, 5), show_value=true)))
	- *distance_decay:* **decay distance** $(Child(:decay, Slider(0.5:0.5:50; default=snap(c.decay_km, 0.5:0.5:50), show_value=true))) km
	- *random:* **link density** $(Child(:density, Slider(0:0.05:1; default=snap(c.density, 0:0.05:1), show_value=true)))

	ADRIA's own model only uses each sink's *relative* mix of sources; absolute values matter for decision criteria and the other engines.
	"""
end

# ╔═╡ 65eedb10-c20a-11f1-3234-255755d44dbf
@bind cover_ui PlutoUI.combine() do Child
	cv = base.cover
	w = cv.composition ./ sum(cv.composition)
	slider(i) = Child(Symbol("w$i"), Slider(0:0.01:1; default=snap(w[i], 0:0.01:1), show_value=true))
	md"""
	## Initial cover
	- **Total cover** (%) $(Child(:total, range_input(cv.total_percent, 0:1:100))) of the $(Child(:basis, Select(["habitable" => "habitable area", "reef" => "whole reef area"]; default=cv.basis)))
	- **Composition** (relative weights)
	  - $(GROUP_LABELS[1]) $(slider(1))
	  - $(GROUP_LABELS[2]) $(slider(2))
	  - $(GROUP_LABELS[3]) $(slider(3))
	  - $(GROUP_LABELS[4]) $(slider(4))
	  - $(GROUP_LABELS[5]) $(slider(5))
	- **Per-reef variation of the mix** $(Child(:jitter, Slider(0:0.01:1; default=snap(cv.jitter, 0:0.01:1), show_value=true))) (0 = same mix everywhere)
	"""
end

# ╔═╡ 65eedb10-c20a-11f1-27a2-2761a4713fc0
@bind dhw_ui PlutoUI.combine() do Child
	d = base.dhw
	labels = sort(collect(keys(d.climate)))
	climate_text = join(["$l = $(SRD.climate_params(d, l).trend_per_decade)" for l in labels], "; ")
	exposure = SRD.attr_spec(d.reef_exposure)
	mode = exposure isa NamedTuple && haskey(exposure, :from) ? "gradient" : "range"
	elo, ehi = mode == "gradient" ? (exposure.from, exposure.to) : SRD.attr_bounds(exposure)
	md"""
	## Heat stress (DHW)
	- **Climate levels** as `label = warming trend (DHW per decade)`, separated by `;` $(Child(:climate, TextField(40; default=climate_text)))
	  Each label writes `DHWs/dhwRCP<label>.nc`. Labels are positive integers that ADRIA treats as numbers, not real RCPs; load with `ADRIA.load_domain(path, label)`. Other per-level overrides in the starting file are kept.
	- **Baseline** $(Child(:baseline, Slider(0:0.1:15; default=snap(d.baseline, 0:0.1:15), show_value=true))) DHW in the first year
	- **Year-to-year variability** SD $(Child(:sd, Slider(0:0.1:6; default=snap(d.interannual_sd, 0:0.1:6), show_value=true))) autocorrelation $(Child(:ar1, Slider(0:0.05:0.95; default=snap(d.ar1, 0:0.05:0.95), show_value=true)))
	- **Marine heatwaves** probability per year $(Child(:event_prob, Slider(0:0.01:1; default=snap(d.event_prob, 0:0.01:1), show_value=true))) adding $(Child(:mag, range_input(d.event_magnitude, 0:0.5:30))) DHW; forced in years $(Child(:event_years, TextField(default=join(d.event_years, ", "))))
	- **Reef exposure** multiplier $(Child(:exp_mode, Select(["range" => "random per reef", "gradient" => "gradient, first → last reef"]; default=mode))) from $(Child(:exp_lo, NumberField(0:0.05:3; default=snap(elo, 0:0.05:3)))) to $(Child(:exp_hi, NumberField(0:0.05:3; default=snap(ehi, 0:0.05:3))))
	- **Reef offset** $(Child(:offset, Slider(-5:0.5:5; default=snap(first(SRD.attr_bounds(d.reef_offset)), -5:0.5:5), show_value=true))) DHW
	- **Local noise** SD $(Child(:local_sd, Slider(0:0.1:6; default=snap(d.local_sd, 0:0.1:6), show_value=true))) (0, with exposure 1 and offset 0, gives identical DHW on every reef)
	"""
end

# ╔═╡ 65eedb10-c20a-11f1-2327-5fc267322216
@bind cyc_ui PlutoUI.combine() do Child
	cy = base.cyclones
	v = cy.group_vulnerability
	slider(i) = Child(Symbol("v$i"), Slider(0:0.05:1; default=snap(v[i], 0:0.05:1), show_value=true))
	md"""
	## Cyclones (system-wide)
	In a cyclone year every reef loses the same fraction of each functional group.
	- $(Child(:enabled, CheckBox(default=cy.enabled))) cyclones enabled
	- **Probability per year** $(Child(:prob, Slider(0:0.01:1; default=snap(cy.annual_prob, 0:0.01:1), show_value=true))); forced in years $(Child(:event_years, TextField(default=join(cy.event_years, ", "))))
	- **Event intensity** (0-1) $(Child(:intensity, range_input(cy.intensity, 0:0.05:1)))
	- **Mortality at intensity 1**
	  - $(GROUP_LABELS[1]) $(slider(1))
	  - $(GROUP_LABELS[2]) $(slider(2))
	  - $(GROUP_LABELS[3]) $(slider(3))
	  - $(GROUP_LABELS[4]) $(slider(4))
	  - $(GROUP_LABELS[5]) $(slider(5))
	"""
end

# ╔═╡ 65eedb10-c20a-11f1-3d5c-5109ea9d4bdc
cfg = try
	let n = reef_ui.n, b = base
	# Per-reef lists from the starting file win over the range controls (if they fit)
	pick(spec, widget) = keep_lists && spec isa AbstractVector && length(spec) == n ? spec : widget
	exposure = if dhw_ui.exp_mode == "gradient"
		(from=Float64(dhw_ui.exp_lo), to=Float64(dhw_ui.exp_hi))
	else
		range_spec((min=dhw_ui.exp_lo, max=dhw_ui.exp_hi))
	end

	DomainConfig(;
		name=String(strip(dom_ui.name)),
		timeframe=[dom_ui.y0, dom_ui.y1],
		seed=dom_ui.seed,
		n_dhw_scenarios=dom_ui.n_dhw,
		n_cyclone_scenarios=dom_ui.n_cyc,
		reefs=ReefConfig(;
			n=n,
			layout=reef_ui.layout,
			center=[reef_ui.lon, reef_ui.lat],
			spacing_km=reef_ui.spacing,
			bearing_deg=reef_ui.bearing,
			n_clusters=min(reef_ui.n_clusters, n),
			cluster_spacing_km=reef_ui.cluster_spacing,
			area_m2=pick(b.reefs.area_m2, range_spec(reef_ui.area)),
			k_percent=pick(b.reefs.k_percent, range_spec(reef_ui.k)),
			depth_med=pick(b.reefs.depth_med, range_spec(reef_ui.depth)),
			ub_med=pick(b.reefs.ub_med, range_spec(reef_ui.ub)),
			habitat=pick(b.reefs.habitat, reef_ui.habitat),
			cb_calib_group=pick(b.reefs.cb_calib_group, reef_ui.cb),
			distinct_ub_depth=b.reefs.distinct_ub_depth,
			waves_from_ub_med=reef_ui.waves_from_ub
		),
		connectivity=ConnectivityConfig(;
			topology=conn_ui.topology,
			self_recruitment=conn_ui.self,
			strength=conn_ui.strength,
			asymmetry=conn_ui.asymmetry,
			reach=conn_ui.reach,
			decay_km=conn_ui.decay,
			density=conn_ui.density
		),
		cover=CoverConfig(;
			total_percent=pick(b.cover.total_percent, range_spec(cover_ui.total)),
			basis=cover_ui.basis,
			composition=[cover_ui.w1, cover_ui.w2, cover_ui.w3, cover_ui.w4, cover_ui.w5],
			jitter=cover_ui.jitter
		),
		dhw=DHWConfig(;
			baseline=dhw_ui.baseline,
			trend_per_decade=b.dhw.trend_per_decade,
			interannual_sd=dhw_ui.sd,
			ar1=dhw_ui.ar1,
			event_prob=dhw_ui.event_prob,
			event_magnitude=[dhw_ui.mag.min, dhw_ui.mag.max],
			event_years=parse_years(dhw_ui.event_years),
			reef_exposure=pick(b.dhw.reef_exposure, exposure),
			reef_offset=pick(b.dhw.reef_offset, Float64(dhw_ui.offset)),
			local_sd=dhw_ui.local_sd,
			climate=parse_climate(dhw_ui.climate, b.dhw.climate)
		),
		cyclones=CycloneConfig(;
			enabled=cyc_ui.enabled,
			annual_prob=cyc_ui.prob,
			event_years=parse_years(cyc_ui.event_years),
			intensity=[cyc_ui.intensity.min, cyc_ui.intensity.max],
			group_vulnerability=[cyc_ui.v1, cyc_ui.v2, cyc_ui.v3, cyc_ui.v4, cyc_ui.v5]
		)
	)
	end
catch err
	err
end

# ╔═╡ 65eedb10-c20a-11f1-1f3e-f33b003ff64a
dom = if cfg isa Exception
	cfg
else
	try
		build_domain(cfg)
	catch err
		err
	end
end

# ╔═╡ 65eedb10-c20a-11f1-2c66-35d8c42dd3d7
md"""
## Preview
"""

# ╔═╡ 65eedb10-c20a-11f1-26fe-d192b6007514
if dom isa Exception
	Markdown.MD(Markdown.Admonition(
		"warning", "The configuration needs fixing",
		[Markdown.Code(sprint(showerror, dom))]
	))
else
	dom
end

# ╔═╡ 65eedb10-c20a-11f1-3151-3d1a1ca58cef
begin
	function plot_map(dom)
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
			title="Reefs and connectivity", subtitle="colour = initial cover, size = area, outline = cluster"
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

	function plot_connectivity(dom)
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

	function plot_cover(dom)
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

	function plot_dhw(dom, label, scen)
		yrs = collect(years(dom))
		D, R = dom.dhw[label], dom.dhw_regional[label]
		exposure = dom.reefs.dhw_exposure
		lo, hi = extrema(exposure)
		fig = Figure(; size=(940, 340))
		ax1 = Axis(fig[1, 1]; title="Label $label, scenario $scen (colour = reef exposure)",
			xlabel="year", ylabel="DHW")
		for r in axes(D, 2)
			lines!(ax1, yrs, D[:, r, scen];
				color=exposure[r], colormap=:thermal, colorrange=(lo, hi + 1e-9), linewidth=1)
		end
		lines!(ax1, yrs, R[:, scen]; color=:black, linewidth=2.5, label="regional signal")
		axislegend(ax1; position=:lt)

		ax2 = Axis(fig[1, 2]; title="All scenarios: mean over reefs (median, 10-90%)", xlabel="year")
		for (k, l) in enumerate(climate_labels(dom))
			m = dropdims(mean(dom.dhw[l]; dims=2); dims=2)
			q = [quantile(m[t, :], p) for t in axes(m, 1), p in (0.1, 0.5, 0.9)]
			c = Makie.wong_colors()[mod1(k, 7)]
			band!(ax2, yrs, q[:, 1], q[:, 3]; color=(c, 0.25))
			lines!(ax2, yrs, q[:, 2]; color=c, linewidth=2, label="label $l")
		end
		axislegend(ax2; position=:lt)
		linkyaxes!(ax1, ax2)
		return fig
	end

	function plot_cyclones(dom, scen)
		yrs = collect(years(dom))
		intensity = dom.cyclone_intensity
		cy = dom.config.cyclones
		fig = Figure(; size=(940, 320))
		n_events = vec(count(>(0), intensity; dims=1))
		ax1 = Axis(fig[1, 1];
			title="Scenario $scen: $(n_events[scen]) cyclone(s); mean $(round(mean(n_events); digits=1)) per scenario",
			xlabel="year", ylabel="intensity", limits=(nothing, (0, 1.05)))
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
end

# ╔═╡ 65eedb10-c20a-11f1-297a-b31697957f7e
dom isa SyntheticDomain ? plot_map(dom) : md"*No preview until the configuration is fixed.*"

# ╔═╡ 65eedb10-c20a-11f1-09ae-91f4864e6f97
dom isa SyntheticDomain ? plot_connectivity(dom) : nothing

# ╔═╡ 65eedb10-c20a-11f1-2579-17de32220f6c
dom isa SyntheticDomain ? plot_cover(dom) : nothing

# ╔═╡ 65eedb10-c20a-11f1-113a-6b7f582a39e3
md"""
### Heat stress
Climate label $(@bind dhw_label Select(dom isa SyntheticDomain ? string.(climate_labels(dom)) : ["-"])) scenario $(@bind dhw_scen Slider(1:(dom isa SyntheticDomain ? dom.config.n_dhw_scenarios : 1); show_value=true))
"""

# ╔═╡ 65eedb10-c20a-11f1-32d5-5fe5c4ff3baf
if dom isa SyntheticDomain
	let labels = climate_labels(dom)
		label = something(tryparse(Int, dhw_label), first(labels))
		label in labels || (label = first(labels))
		plot_dhw(dom, label, min(dhw_scen, dom.config.n_dhw_scenarios))
	end
end

# ╔═╡ 65eedb10-c20a-11f1-2c11-d78ab4f982a8
md"""
### Cyclones
Scenario $(@bind cyc_scen Slider(1:(dom isa SyntheticDomain ? dom.config.n_cyclone_scenarios : 1); show_value=true))
"""

# ╔═╡ 65eedb10-c20a-11f1-1870-a3ea1e1ad31a
dom isa SyntheticDomain ? plot_cyclones(dom, min(cyc_scen, dom.config.n_cyclone_scenarios)) : nothing

# ╔═╡ 65eedb10-c20a-11f1-17b4-df5ca3a3fd27
md"""
## Write
- Domain folder $(@bind outdir TextField(70; default=joinpath(REPO, "domains", strip(dom_ui.name))))
- Config file $(@bind config_out TextField(70; default=joinpath(CONFIG_DIR, strip(dom_ui.name) * ".toml")))
- $(@bind overwrite CheckBox(default=false)) replace an existing domain folder / config file

$(@bind write_clicks CounterButton("Write domain")) $(@bind save_clicks CounterButton("Save config"))
"""

# ╔═╡ 65eedb10-c20a-11f1-0024-cdc4f1737d36
# Latest values, held in Refs so the button cells below re-run only when clicked
# (not on every slider move).
latest = Ref{Any}(nothing)

# ╔═╡ 65eedb10-c20a-11f1-1a04-410ba66d28ca
latest[] = (dom=dom, cfg=cfg, outdir=String(strip(outdir)), config_out=String(strip(config_out)), overwrite=overwrite);

# ╔═╡ 65eedb10-c20a-11f1-1e5f-edcd3f50528c
if write_clicks > 0 && latest[] !== nothing
	let s = latest[]
		if s.dom isa SyntheticDomain
			try
				path = write_domain(s.dom, s.outdir; force=s.overwrite)
				label = first(climate_labels(s.dom))
				Markdown.parse("""
				Domain written to `$(path)`. Load it with

				```julia
				dom = ADRIA.load_domain(raw"$(path)", "$(label)")
				```
				""")
			catch err
				Markdown.MD(Markdown.Admonition("danger", "Not written", [Markdown.Paragraph([sprint(showerror, err)])]))
			end
		else
			md"**Not written:** fix the configuration first."
		end
	end
end

# ╔═╡ 65efec80-c20a-11f1-06b7-e1d4016d982e
if save_clicks > 0 && latest[] !== nothing
	let s = latest[]
		if !(s.cfg isa DomainConfig)
			md"**Not saved:** fix the configuration first."
		elseif isfile(s.config_out) && !s.overwrite
			md"**Not saved:** `$(s.config_out)` exists. Tick *replace* to overwrite it."
		else
			save_config(s.config_out, s.cfg)
			Markdown.parse("""
			Config saved to `$(s.config_out)`. Regenerate the same domain with

			```julia
			generate_domain(raw"$(s.config_out)")
			```
			""")
		end
	end
end

# ╔═╡ Cell order:
# ╟─65ed2d60-c20a-11f1-203a-a7b44f43174d
# ╟─65eedb10-c20a-11f1-049d-5511fec48719
# ╟─65eedb10-c20a-11f1-0155-a1518672437e
# ╟─65eedb10-c20a-11f1-0fff-792e6ec4859a
# ╟─65eedb10-c20a-11f1-37f7-652c033e5226
# ╟─65eedb10-c20a-11f1-3a88-459f1b87407e
# ╟─65eedb10-c20a-11f1-37ed-37e36fc02c3c
# ╟─65eedb10-c20a-11f1-3921-d9d49b19e4ae
# ╟─65eedb10-c20a-11f1-3f17-1d1b42bb3704
# ╟─65eedb10-c20a-11f1-3234-255755d44dbf
# ╟─65eedb10-c20a-11f1-27a2-2761a4713fc0
# ╟─65eedb10-c20a-11f1-2327-5fc267322216
# ╟─65eedb10-c20a-11f1-3d5c-5109ea9d4bdc
# ╟─65eedb10-c20a-11f1-1f3e-f33b003ff64a
# ╟─65eedb10-c20a-11f1-2c66-35d8c42dd3d7
# ╟─65eedb10-c20a-11f1-26fe-d192b6007514
# ╟─65eedb10-c20a-11f1-3151-3d1a1ca58cef
# ╟─65eedb10-c20a-11f1-297a-b31697957f7e
# ╟─65eedb10-c20a-11f1-09ae-91f4864e6f97
# ╟─65eedb10-c20a-11f1-2579-17de32220f6c
# ╟─65eedb10-c20a-11f1-113a-6b7f582a39e3
# ╟─65eedb10-c20a-11f1-32d5-5fe5c4ff3baf
# ╟─65eedb10-c20a-11f1-2c11-d78ab4f982a8
# ╟─65eedb10-c20a-11f1-1870-a3ea1e1ad31a
# ╟─65eedb10-c20a-11f1-17b4-df5ca3a3fd27
# ╟─65eedb10-c20a-11f1-0024-cdc4f1737d36
# ╟─65eedb10-c20a-11f1-1a04-410ba66d28ca
# ╟─65eedb10-c20a-11f1-1e5f-edcd3f50528c
# ╟─65efec80-c20a-11f1-06b7-e1d4016d982e
