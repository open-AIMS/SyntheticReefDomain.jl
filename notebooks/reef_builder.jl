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

# ╔═╡ 30391b60-c2c9-11f1-357b-e56922d03858
begin
	import Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	using SyntheticReefDomain, CairoMakie, Statistics
	# Explicit imports: Makie also exports a `Slider`
	using PlutoUI: PlutoUI, Slider, Select, CheckBox, TextField, NumberField, CounterButton,
		FilePicker, DownloadButton
	import SyntheticReefDomain as SRD
	CairoMakie.activate!(; type="png", px_per_unit=2)
	const REPO = normpath(joinpath(@__DIR__, ".."))
	const CONFIG_DIR = joinpath(REPO, "configs")
end

# ╔═╡ 3036d170-c2c9-11f1-3da2-310dc40f6538
md"""
# Synthetic reef builder

Build a synthetic ADRIA domain with 1-20 reefs. Pick a starting configuration, adjust the
reef system, connectivity, initial cover and disturbances, check the previews, then
**write the domain** -- or save the settings as a TOML config and regenerate exactly the
same domain from a script with `generate_domain("configs/<name>.toml")`.

Every change rebuilds the domain in memory; nothing is written to disk until you press one
of the buttons at the bottom.
"""

# ╔═╡ 30391b60-c2c9-11f1-26e7-314bdf737a3f
config_files = let
	# Configs in configs/ and configs/examples/, plus the settings of every domain already
	# written to domains/ (each domain folder holds its generator_config.toml)
	files = String[]
	for dir in (CONFIG_DIR, joinpath(CONFIG_DIR, "examples"))
		isdir(dir) && append!(files, [joinpath(dir, f) for f in readdir(dir) if endswith(f, ".toml")])
	end
	domains_dir = joinpath(REPO, "domains")
	if isdir(domains_dir)
		for d in readdir(domains_dir; join=true)
			f = joinpath(d, "generator_config.toml")
			isfile(f) && push!(files, f)
		end
	end
	label(f) = replace(relpath(f, REPO), '\\' => '/')
	sort!(files; by=f -> (label(f) != "configs/default.toml", !startswith(label(f), "configs/"), label(f)))
	[f => label(f) for f in files]
end

# ╔═╡ 30391b60-c2c9-11f1-13e3-8d2c4c1e5b6e
md"""
## Starting point
Start from a config in this repository (`configs/`, or a domain already written to `domains/`): $(@bind base_file Select(config_files))
"""

# ╔═╡ 30391b60-c2c9-11f1-2a7c-f910d1de1b83
begin
	base_file  # choosing from the list above clears a file opened below
	md"""
	... or open any config file (opens the file dialog): $(@bind picked_file FilePicker([MIME(".toml")]))
	"""
end

# ╔═╡ 30391b60-c2c9-11f1-0556-05b886114c8a
begin
	# A file opened with the picker wins over the list; an unreadable file falls back to it
	list_choice() = (cfg=load_config(base_file),
		source="`$(replace(relpath(base_file, REPO), '\\' => '/'))`", error=nothing)
	base_choice = if picked_file === nothing
		list_choice()
	else
		try
			(cfg=load_config(IOBuffer(picked_file["data"])),
				source="`$(picked_file["name"])` (opened file)", error=nothing)
		catch err
			merge(list_choice(), (error=sprint(showerror, err),))
		end
	end
	base = base_choice.cfg
end

# ╔═╡ 30391b60-c2c9-11f1-3290-f1c758b41acd
if base_choice.error === nothing
	Markdown.parse("Starting values from $(base_choice.source).")
else
	Markdown.MD(Markdown.Admonition("warning", "Could not use the opened file", [
		Markdown.Paragraph(["It is not a valid configuration, so the starting values come from " *
			"$(base_choice.source) instead:"]),
		Markdown.Code(base_choice.error)
	]))
end

# ╔═╡ 30391b60-c2c9-11f1-122c-b14f46a9fee0
md"""
$(@bind keep_lists CheckBox(default=true)) keep per-reef lists from the starting file (used
when their length matches the number of reefs; the matching range controls below are then
ignored).
"""

# ╔═╡ 30391b60-c2c9-11f1-0974-cbe7b3d5cbee
begin
	"Closest value of `r` to `v` (widget defaults must lie on the widget's grid)."
	snap(v, r) = r[argmin(abs.(r .- v))]

	# Min and max number fields initialised from a per-reef attribute specification.
	# (Kept as two separate children: a `combine` nested inside another `combine` sends
	# `nothing` from the browser.)
	min_field(spec, r) = NumberField(r; default=snap(first(SRD.attr_bounds(spec)), r))
	max_field(spec, r) = NumberField(r; default=snap(last(SRD.attr_bounds(spec)), r))

	"Per-reef specification from a min and a max (a single value when they are equal)."
	function range_spec(a, b)
		lo, hi = minmax(Float64(a), Float64(b))
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
end

# ╔═╡ 30391b60-c2c9-11f1-0e5e-8160066b836c
@bind dom_ui PlutoUI.combine() do Child
	md"""
	## Domain
	- **Name** (also the folder name) $(Child(:name, TextField(default=base.name)))
	- **Years** $(Child(:y0, NumberField(1950:2300; default=base.timeframe[1]))) to $(Child(:y1, NumberField(1950:2300; default=base.timeframe[2]))) (at least 26 years)
	- **Seed** $(Child(:seed, NumberField(0:1_000_000_000; default=base.seed)))
	- Stochastic scenarios (every file also gets scenario 1 with no disturbance):
	- **DHW scenarios** $(Child(:n_dhw, Slider(1:100; default=clamp(base.n_dhw_scenarios, 1, 100), show_value=true)))
	- **Cyclone scenarios** $(Child(:n_cyc, Slider(1:100; default=clamp(base.n_cyclone_scenarios, 1, 100), show_value=true)))
	"""
end

# ╔═╡ 30391b60-c2c9-11f1-10ca-61570837d1a9
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
	- **Area** (m²) $(Child(:area_min, min_field(r.area_m2, 1000:1000:5_000_000))) to $(Child(:area_max, max_field(r.area_m2, 1000:1000:5_000_000)))
	- **Habitable fraction k** (%) $(Child(:k_min, min_field(r.k_percent, 1:100))) to $(Child(:k_max, max_field(r.k_percent, 1:100)))
	- **Median depth** (m) $(Child(:depth_min, min_field(r.depth_med, 0.5:0.5:40))) to $(Child(:depth_max, max_field(r.depth_med, 0.5:0.5:40)))
	- **Median wave orbital velocity `ub_med`** (m/s, read by C~scape) $(Child(:ub_min, min_field(r.ub_med, 0:0.01:2))) to $(Child(:ub_max, max_field(r.ub_med, 0:0.01:2)))
	- **Habitat** $(Child(:habitat, Select(HABITATS; default=habitat))) **CoralBlox calibration group** $(Child(:cb, NumberField(1:12; default=cb)))
	- $(Child(:waves_from_ub, CheckBox(default=r.waves_from_ub_med))) write each reef's `ub_med` as its (constant) wave series
	"""
end

# ╔═╡ 30391b60-c2c9-11f1-3ad1-d7ce51a8e5d9
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
	- *truncated_dilution* (distance-based connectivity for C~scape; self-recruitment and strength are not used): **no settlement within** $(Child(:threshold, NumberField(0:0.01:2000; default=snap(c.threshold_m, 0:0.01:2000)))) m of a source reef's centre, **dispersal distance** (kernel down to 5%) $(Child(:dispersal, Slider(1:100; default=snap(c.dispersal_km, 1:100), show_value=true))) km, **matrix total** $(Child(:matrix_total, NumberField(0.00001:0.00001:1; default=snap(c.matrix_total, 0.00001:0.00001:1))))

	ADRIA's own model only uses each sink's *relative* mix of sources; absolute values matter for decision criteria and the other engines.
	"""
end

# ╔═╡ 30391b60-c2c9-11f1-0de7-79e7d624caac
@bind cover_ui PlutoUI.combine() do Child
	cv = base.cover
	w = cv.composition ./ sum(cv.composition)
	slider(i) = Child(Symbol("w$i"), Slider(0:0.01:1; default=snap(w[i], 0:0.01:1), show_value=true))
	md"""
	## Initial cover
	- **Total cover** (%) $(Child(:total_min, min_field(cv.total_percent, 0:1:100))) to $(Child(:total_max, max_field(cv.total_percent, 0:1:100))) of the $(Child(:basis, Select(["habitable" => "habitable area", "reef" => "whole reef area"]; default=cv.basis)))
	- **Composition** (relative weights)
	  - $(GROUP_LABELS[1]) $(slider(1))
	  - $(GROUP_LABELS[2]) $(slider(2))
	  - $(GROUP_LABELS[3]) $(slider(3))
	  - $(GROUP_LABELS[4]) $(slider(4))
	  - $(GROUP_LABELS[5]) $(slider(5))
	- **Per-reef variation of the mix** $(Child(:jitter, Slider(0:0.01:1; default=snap(cv.jitter, 0:0.01:1), show_value=true))) (0 = same mix everywhere)
	"""
end

# ╔═╡ 30391b60-c2c9-11f1-2f55-bda5e807038f
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
	- **Marine heatwaves** probability per year $(Child(:event_prob, Slider(0:0.01:1; default=snap(d.event_prob, 0:0.01:1), show_value=true))) adding $(Child(:mag_min, min_field(d.event_magnitude, 0:0.5:30))) to $(Child(:mag_max, max_field(d.event_magnitude, 0:0.5:30))) DHW; forced in years $(Child(:event_years, TextField(default=join(d.event_years, ", "))))
	- **Reef exposure** multiplier $(Child(:exp_mode, Select(["range" => "random per reef", "gradient" => "gradient, first → last reef"]; default=mode))) from $(Child(:exp_lo, NumberField(0:0.05:3; default=snap(elo, 0:0.05:3)))) to $(Child(:exp_hi, NumberField(0:0.05:3; default=snap(ehi, 0:0.05:3))))
	- **Reef offset** $(Child(:offset, Slider(-5:0.5:5; default=snap(first(SRD.attr_bounds(d.reef_offset)), -5:0.5:5), show_value=true))) DHW
	- **Local noise** SD $(Child(:local_sd, Slider(0:0.1:6; default=snap(d.local_sd, 0:0.1:6), show_value=true))) (0, with exposure 1 and offset 0, gives identical DHW on every reef)
	"""
end

# ╔═╡ 30391b60-c2c9-11f1-1066-d3b3fbfd9aa7
@bind cyc_ui PlutoUI.combine() do Child
	cy = base.cyclones
	v = cy.group_vulnerability
	slider(i) = Child(Symbol("v$i"), Slider(0:0.05:1; default=snap(v[i], 0:0.05:1), show_value=true))
	md"""
	## Cyclones (system-wide)
	In a cyclone year every reef loses the same fraction of each functional group.
	- $(Child(:enabled, CheckBox(default=cy.enabled))) cyclones enabled
	- **Probability per year** $(Child(:prob, Slider(0:0.01:1; default=snap(cy.annual_prob, 0:0.01:1), show_value=true))); forced in years $(Child(:event_years, TextField(default=join(cy.event_years, ", "))))
	- **Event intensity** (0-1) $(Child(:int_min, min_field(cy.intensity, 0:0.05:1))) to $(Child(:int_max, max_field(cy.intensity, 0:0.05:1)))
	- **Mortality at intensity 1**
	  - $(GROUP_LABELS[1]) $(slider(1))
	  - $(GROUP_LABELS[2]) $(slider(2))
	  - $(GROUP_LABELS[3]) $(slider(3))
	  - $(GROUP_LABELS[4]) $(slider(4))
	  - $(GROUP_LABELS[5]) $(slider(5))
	"""
end

# ╔═╡ 30391b60-c2c9-11f1-058a-6f7e36312428
cfg = try
	let n = reef_ui.n, b = base
	# Per-reef lists from the starting file win over the range controls (if they fit)
	pick(spec, widget) = keep_lists && spec isa AbstractVector && length(spec) == n ? spec : widget
	exposure = if dhw_ui.exp_mode == "gradient"
		(from=Float64(dhw_ui.exp_lo), to=Float64(dhw_ui.exp_hi))
	else
		range_spec(dhw_ui.exp_lo, dhw_ui.exp_hi)
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
			area_m2=pick(b.reefs.area_m2, range_spec(reef_ui.area_min, reef_ui.area_max)),
			k_percent=pick(b.reefs.k_percent, range_spec(reef_ui.k_min, reef_ui.k_max)),
			depth_med=pick(b.reefs.depth_med, range_spec(reef_ui.depth_min, reef_ui.depth_max)),
			ub_med=pick(b.reefs.ub_med, range_spec(reef_ui.ub_min, reef_ui.ub_max)),
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
			density=conn_ui.density,
			threshold_m=Float64(conn_ui.threshold),
			dispersal_km=Float64(conn_ui.dispersal),
			matrix_total=Float64(conn_ui.matrix_total)
		),
		cover=CoverConfig(;
			total_percent=pick(b.cover.total_percent, range_spec(cover_ui.total_min, cover_ui.total_max)),
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
			event_magnitude=sort([dhw_ui.mag_min, dhw_ui.mag_max]),
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
			intensity=sort([cyc_ui.int_min, cyc_ui.int_max]),
			group_vulnerability=[cyc_ui.v1, cyc_ui.v2, cyc_ui.v3, cyc_ui.v4, cyc_ui.v5]
		)
	)
	end
catch err
	err
end

# ╔═╡ 30391b60-c2c9-11f1-0e8f-7df8023603df
dom = if cfg isa Exception
	cfg
else
	try
		build_domain(cfg)
	catch err
		err
	end
end

# ╔═╡ 30391b60-c2c9-11f1-3492-c575ff1eff9f
md"""
## Preview
"""

# ╔═╡ 30391b60-c2c9-11f1-2887-694d8bc7e804
if dom isa Exception
	Markdown.MD(Markdown.Admonition(
		"warning", "The configuration needs fixing",
		[Markdown.Code(sprint(showerror, dom))]
	))
else
	dom
end

# ╔═╡ 30391b60-c2c9-11f1-2680-d5e27861cb41
dom isa SyntheticDomain ? plot_map(dom) : md"*No preview until the configuration is fixed.*"

# ╔═╡ 30391b60-c2c9-11f1-2437-69d5671aabbb
dom isa SyntheticDomain ? plot_connectivity(dom) : nothing

# ╔═╡ 30391b60-c2c9-11f1-15df-87704e7828a2
dom isa SyntheticDomain ? plot_cover(dom) : nothing

# ╔═╡ 30391b60-c2c9-11f1-08af-9d179115b5a8
md"""
### Heat stress
Climate label $(@bind dhw_label Select(dom isa SyntheticDomain ? string.(climate_labels(dom)) : ["-"])) scenario (1 = no heat stress) $(@bind dhw_scen Slider(1:(dom isa SyntheticDomain ? size(first(values(dom.dhw)), 3) : 2); default=2, show_value=true))
"""

# ╔═╡ 30391b60-c2c9-11f1-3ba5-7b9850f61436
if dom isa SyntheticDomain
	let labels = climate_labels(dom)
		label = something(tryparse(Int, dhw_label), first(labels))
		label in labels || (label = first(labels))
		plot_dhw(dom, label, min(dhw_scen, size(dom.dhw[label], 3)))
	end
end

# ╔═╡ 30391b60-c2c9-11f1-31b6-7984b6b2a5e3
md"""
### Cyclones
Scenario (1 = no cyclones) $(@bind cyc_scen Slider(1:(dom isa SyntheticDomain ? size(dom.cyclone_intensity, 2) : 2); default=2, show_value=true))
"""

# ╔═╡ 30391b60-c2c9-11f1-01a6-679afafb69b6
dom isa SyntheticDomain ? plot_cyclones(dom, min(cyc_scen, size(dom.cyclone_intensity, 2))) : nothing

# ╔═╡ 30391b60-c2c9-11f1-1538-41d5b5f7e277
md"""
## Write
- Domain folder $(@bind outdir TextField(70; default=joinpath(REPO, "domains", strip(dom_ui.name))))
- Config file $(@bind config_out TextField(70; default=joinpath(CONFIG_DIR, strip(dom_ui.name) * ".toml")))
- $(@bind overwrite CheckBox(default=false)) replace an existing domain folder / config file

$(@bind write_clicks CounterButton("Write domain")) $(@bind save_clicks CounterButton("Save config")) $(@bind zip_clicks CounterButton("Prepare download (.zip)"))
"""

# ╔═╡ 30391b60-c2c9-11f1-0a45-dd1365e29010
# Latest values, held in Refs so the button cells below re-run only when clicked
# (not on every slider move).
latest = Ref{Any}(nothing)

# ╔═╡ 30391b60-c2c9-11f1-1887-b3ca9a7b5fa3
latest[] = (dom=dom, cfg=cfg, outdir=String(strip(outdir)), config_out=String(strip(config_out)), overwrite=overwrite);

# ╔═╡ 30391b60-c2c9-11f1-34d9-67496ba3ff08
if write_clicks > 0 && latest[] !== nothing
	let s = latest[]
		if s.dom isa SyntheticDomain
			try
				path = write_domain(s.dom, s.outdir; force=s.overwrite)
				label = first(climate_labels(s.dom))
				Markdown.parse("""
				Domain written to `$(path)`. Its `README.md` lists every setting and shows these plots (saved in `figures/`). Load it with

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

# ╔═╡ 30391b60-c2c9-11f1-33ab-b9bd2a11a812
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

# ╔═╡ 30391b60-c2c9-11f1-057f-b1ab42c77c5a
if zip_clicks > 0 && latest[] !== nothing
	let s = latest[]
		if s.dom isa SyntheticDomain
			DownloadButton(domain_zip(s.dom), "$(s.dom.config.name).zip")
		else
			md"**No download:** fix the configuration first."
		end
	end
end

# ╔═╡ Cell order:
# ╟─3036d170-c2c9-11f1-3da2-310dc40f6538
# ╟─30391b60-c2c9-11f1-357b-e56922d03858
# ╟─30391b60-c2c9-11f1-26e7-314bdf737a3f
# ╟─30391b60-c2c9-11f1-13e3-8d2c4c1e5b6e
# ╟─30391b60-c2c9-11f1-2a7c-f910d1de1b83
# ╟─30391b60-c2c9-11f1-0556-05b886114c8a
# ╟─30391b60-c2c9-11f1-3290-f1c758b41acd
# ╟─30391b60-c2c9-11f1-122c-b14f46a9fee0
# ╟─30391b60-c2c9-11f1-0974-cbe7b3d5cbee
# ╟─30391b60-c2c9-11f1-0e5e-8160066b836c
# ╟─30391b60-c2c9-11f1-10ca-61570837d1a9
# ╟─30391b60-c2c9-11f1-3ad1-d7ce51a8e5d9
# ╟─30391b60-c2c9-11f1-0de7-79e7d624caac
# ╟─30391b60-c2c9-11f1-2f55-bda5e807038f
# ╟─30391b60-c2c9-11f1-1066-d3b3fbfd9aa7
# ╟─30391b60-c2c9-11f1-058a-6f7e36312428
# ╟─30391b60-c2c9-11f1-0e8f-7df8023603df
# ╟─30391b60-c2c9-11f1-3492-c575ff1eff9f
# ╟─30391b60-c2c9-11f1-2887-694d8bc7e804
# ╟─30391b60-c2c9-11f1-2680-d5e27861cb41
# ╟─30391b60-c2c9-11f1-2437-69d5671aabbb
# ╟─30391b60-c2c9-11f1-15df-87704e7828a2
# ╟─30391b60-c2c9-11f1-08af-9d179115b5a8
# ╟─30391b60-c2c9-11f1-3ba5-7b9850f61436
# ╟─30391b60-c2c9-11f1-31b6-7984b6b2a5e3
# ╟─30391b60-c2c9-11f1-01a6-679afafb69b6
# ╟─30391b60-c2c9-11f1-1538-41d5b5f7e277
# ╟─30391b60-c2c9-11f1-0a45-dd1365e29010
# ╟─30391b60-c2c9-11f1-1887-b3ca9a7b5fa3
# ╟─30391b60-c2c9-11f1-34d9-67496ba3ff08
# ╟─30391b60-c2c9-11f1-33ab-b9bd2a11a812
# ╟─30391b60-c2c9-11f1-057f-b1ab42c77c5a
