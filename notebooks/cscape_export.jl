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

# ╔═╡ 32de4ac0-c2c9-11f1-3937-bb5adc4f7077
begin
	import Pkg
	Pkg.activate(@__DIR__)
	Pkg.instantiate()
	using SyntheticReefDomain, CairoMakie
	# Explicit imports: Makie also exports some of these names
	using PlutoUI: PlutoUI, Select, CheckBox, TextField, NumberField, CounterButton,
		FilePicker, DownloadButton
	import SyntheticReefDomain as SRD
	CairoMakie.activate!(; type="png", px_per_unit=2)
	const REPO = normpath(joinpath(@__DIR__, ".."))
	const CONFIG_DIR = joinpath(REPO, "configs")
end

# ╔═╡ 32db1670-c2c9-11f1-26b9-7d8e198b2e41
md"""
# C~scape (R) export

Turn a synthetic domain into an input project for the R package
[cscape](https://github.com/open-AIMS/c-scape), which runs with `run_cscape()`. Pick a domain:
its settings are rebuilt here, so the C~scape inputs match the ADRIA domain exactly. Then
set the C~scape options, check the previews, and write the project or download it as a zip.

Demography files (`demog_inputs_<type>.rds`) are **not** created: add them to the project's
`data/` folder before running C~scape.
"""

# ╔═╡ 32de4ac0-c2c9-11f1-2cc9-d9d7308631e9
config_files = let
	# Domains already written to domains/ first, then the configs in configs/
	files = String[]
	domains_dir = joinpath(REPO, "domains")
	if isdir(domains_dir)
		for d in readdir(domains_dir; join=true)
			f = joinpath(d, "generator_config.toml")
			isfile(f) && push!(files, f)
		end
	end
	for dir in (CONFIG_DIR, joinpath(CONFIG_DIR, "examples"))
		isdir(dir) && append!(files, sort([joinpath(dir, f) for f in readdir(dir) if endswith(f, ".toml")]))
	end
	label(f) = replace(relpath(f, REPO), '\\' => '/')
	[f => label(f) for f in files]
end

# ╔═╡ 32de4ac0-c2c9-11f1-3a46-fba2ed80b045
md"""
## Domain
Domain to export (written to `domains/`, or a config in `configs/`): $(@bind base_file Select(config_files))
"""

# ╔═╡ 32de4ac0-c2c9-11f1-3a54-cfe7b33e29a8
begin
	base_file  # choosing from the list above clears a file opened below
	md"""
	... or open any config file, e.g. a domain's `generator_config.toml` (opens the file dialog): $(@bind picked_file FilePicker([MIME(".toml")]))
	"""
end

# ╔═╡ 32de4ac0-c2c9-11f1-0de8-83b248e41252
begin
	label_of(f) = replace(relpath(f, REPO), '\\' => '/')
	from_list() = (cfg=load_config(base_file), source="`$(label_of(base_file))`",
		domain_dir=basename(base_file) == "generator_config.toml" ? dirname(base_file) : nothing,
		error=nothing)
	base_choice = if picked_file === nothing
		from_list()
	else
		try
			(cfg=load_config(IOBuffer(picked_file["data"])),
				source="`$(picked_file["name"])` (opened file)", domain_dir=nothing, error=nothing)
		catch err
			merge(from_list(), (error=sprint(showerror, err),))
		end
	end
	base = base_choice.cfg
	# Write into the domain's own folder when exporting a domain from domains/
	default_target = joinpath(something(base_choice.domain_dir, joinpath(REPO, "domains", base.name)), "cscape")
end

# ╔═╡ 32de4ac0-c2c9-11f1-32c2-074f4ceac6fa
if base_choice.error === nothing
	Markdown.parse("Exporting domain **$(base.name)** from $(base_choice.source).")
else
	Markdown.MD(Markdown.Admonition("warning", "Could not use the opened file", [
		Markdown.Paragraph(["It is not a valid configuration, so $(base_choice.source) is used instead:"]),
		Markdown.Code(base_choice.error)
	]))
end

# ╔═╡ 32de4ac0-c2c9-11f1-0f09-4b801ff9b486
dom = try
	build_domain(base)
catch err
	err
end

# ╔═╡ 32de4ac0-c2c9-11f1-32fd-9942596e0877
if dom isa Exception
	Markdown.MD(Markdown.Admonition("warning", "The domain configuration needs fixing",
		[Markdown.Code(sprint(showerror, dom))]))
else
	dom
end

# ╔═╡ 32de4ac0-c2c9-11f1-2640-658b14be5b60
# Start from the options of an earlier export of this domain, if there is one
start_options = let f = joinpath(default_target, "cscape_export.toml")
	isfile(f) ? (try load_cscape_options(f) catch; CscapeOptions() end) : CscapeOptions()
end

# ╔═╡ 32de4ac0-c2c9-11f1-2e56-4b36878aa0d4
@bind opt_ui PlutoUI.combine() do Child
	o = start_options
	region = o.region in CSCAPE_REGIONS ? o.region : first(CSCAPE_REGIONS)
	md"""
	## C~scape options
	- **Region** (a `Cluster` level of the demography models) $(Child(:region, Select(CSCAPE_REGIONS; default=region)))
	- **Demography file names**, one per line in this order: $(join(CSCAPE_TYPES, ", ")) $(Child(:demog, TextField((60, 6); default=join(o.demography_files, "\n"))))
	- **Spawner connectivity scale** $(Child(:scale, NumberField(0:0.0001:10; default=o.connectivity_scale))) (C~scape uses absolute transfer probabilities; ADRIA only uses proportions)
	- **Brooder**: dispersal per source reef $(Child(:brood, NumberField(0:0.001:1; default=o.brood_site_value))) within $(Child(:radius, NumberField(0:0.1:50; default=o.brooder_neighbour_km))) km (90% stays home)
	"""
end

# ╔═╡ 32de4ac0-c2c9-11f1-29e1-4f6aae6f8d36
opts = try
	SRD.validate(CscapeOptions(;
		region=opt_ui.region,
		demography_files=String[strip(l) for l in split(opt_ui.demog, '\n') if !isempty(strip(l))],
		connectivity_scale=Float64(opt_ui.scale),
		brood_site_value=Float64(opt_ui.brood),
		brooder_neighbour_km=Float64(opt_ui.radius)
	))
catch err
	err
end

# ╔═╡ 32de4ac0-c2c9-11f1-2151-a1386d958103
inputs = if dom isa Exception || opts isa Exception
	nothing
else
	try
		cscape_inputs(dom, opts)
	catch err
		err
	end
end

# ╔═╡ 32de4ac0-c2c9-11f1-132c-811be0ca02c5
let
	status = rscript_status()
	r_line = if status.rscript === nothing
		"**R not found:** the export writes the CSV inputs and `build_cscape_project.R`. Run that script later on a machine with R and cscape to create the `.rds` files."
	elseif !status.cscape
		"**R found, but the cscape package is not installed:** the export writes the inputs and the build script. Install cscape, then run the script."
	else
		"**R and cscape found:** the export builds the `.rds` files and `ScenarioID.xlsx` straight away."
	end
	blocks = Any[Markdown.parse(r_line)]
	opts isa Exception && push!(blocks, Markdown.Admonition("warning", "Options need fixing", [Markdown.Code(sprint(showerror, opts))]))
	inputs isa Exception && push!(blocks, Markdown.Admonition("warning", "Could not translate the domain", [Markdown.Code(sprint(showerror, inputs))]))
	if inputs isa NamedTuple && !isempty(inputs.warnings)
		push!(blocks, Markdown.Admonition("note", "Outside the demography models' range",
			[Markdown.Paragraph([w]) for w in inputs.warnings]))
	end
	push!(blocks, Markdown.parse("""
	How the domain is translated: per-reef DHW goes to `DHW_site`; cyclones become categories 1-5 (`ceil(5 x intensity)`), the same on every reef, and C~scape applies its own wind-speed mortality (the domain's per-group cyclone vulnerability is not used); replicate *r* uses DHW scenario *r*, so **replicate 1 is the disturbance-free control**; no COTS and no external eggs.
	"""))
	Markdown.MD(blocks)
end

# ╔═╡ 32de4ac0-c2c9-11f1-1109-413fcf0723ee
md"""
## Preview
Site table as written for C~scape (`k` in %, `area` in m², initial cover in % of site area):
"""

# ╔═╡ 32de4ac0-c2c9-11f1-1028-e39ac09ff156
inputs isa NamedTuple ? hcat(inputs.sites, inputs.init_cover[:, 2:end]; makeunique=true) : nothing

# ╔═╡ 32de4ac0-c2c9-11f1-3acb-2f01b10eb837
inputs isa NamedTuple ? plot_cscape_connectivity(dom, opts) : nothing

# ╔═╡ 32de4ac0-c2c9-11f1-0f14-95f84c34976c
inputs isa NamedTuple ? plot_cscape_cyclones(dom, opts) : nothing

# ╔═╡ 32de4ac0-c2c9-11f1-3e7b-5d2fca3bf741
inputs isa NamedTuple ? plot_cscape_cover(dom, opts) : nothing

# ╔═╡ 32de4ac0-c2c9-11f1-1dc7-fbb2043060d6
md"""
## Write
- Project folder $(@bind target TextField(70; default=default_target))
- $(@bind overwrite CheckBox(default=false)) update an existing export (files you added there, such as demography files, are kept)

$(@bind write_clicks CounterButton("Write C~scape project")) $(@bind zip_clicks CounterButton("Prepare download (.zip)"))
"""

# ╔═╡ 32de4ac0-c2c9-11f1-22a3-89e1202105ea
# Latest values, held in a Ref so the button cells below re-run only when clicked
latest = Ref{Any}(nothing)

# ╔═╡ 32de4ac0-c2c9-11f1-3910-715c854fa3ba
latest[] = (dom=dom, opts=opts, inputs=inputs, target=String(strip(target)), overwrite=overwrite);

# ╔═╡ 32de4ac0-c2c9-11f1-1eb9-65b38f857339
if write_clicks > 0 && latest[] !== nothing
	let s = latest[]
		if !(s.inputs isa NamedTuple)
			md"**Not written:** fix the problems above first."
		else
			try
				r = write_cscape_project(s.dom, s.target; options=s.opts, force=s.overwrite)
				label = first(climate_labels(s.dom))
				fpath = replace(r.path, '\\' => '/')
				Markdown.parse("""
				**$(r.built ? "Built" : "Written"):** `$(r.path)`. $(r.message)

				Add the demography files to `data/`, then run in R:

				```r
				library(cscape)
				run_cscape_batch($label, "$fpath")
				```
				""")
			catch err
				Markdown.MD(Markdown.Admonition("danger", "Not written", [Markdown.Paragraph([sprint(showerror, err)])]))
			end
		end
	end
end

# ╔═╡ 32de4ac0-c2c9-11f1-35a9-8f5fef32f61e
if zip_clicks > 0 && latest[] !== nothing
	let s = latest[]
		if s.inputs isa NamedTuple
			DownloadButton(cscape_zip(s.dom, s.opts), "$(s.dom.config.name)_cscape.zip")
		else
			md"**No download:** fix the problems above first."
		end
	end
end

# ╔═╡ Cell order:
# ╟─32db1670-c2c9-11f1-26b9-7d8e198b2e41
# ╟─32de4ac0-c2c9-11f1-3937-bb5adc4f7077
# ╟─32de4ac0-c2c9-11f1-2cc9-d9d7308631e9
# ╟─32de4ac0-c2c9-11f1-3a46-fba2ed80b045
# ╟─32de4ac0-c2c9-11f1-3a54-cfe7b33e29a8
# ╟─32de4ac0-c2c9-11f1-0de8-83b248e41252
# ╟─32de4ac0-c2c9-11f1-32c2-074f4ceac6fa
# ╟─32de4ac0-c2c9-11f1-0f09-4b801ff9b486
# ╟─32de4ac0-c2c9-11f1-32fd-9942596e0877
# ╟─32de4ac0-c2c9-11f1-2640-658b14be5b60
# ╟─32de4ac0-c2c9-11f1-2e56-4b36878aa0d4
# ╟─32de4ac0-c2c9-11f1-29e1-4f6aae6f8d36
# ╟─32de4ac0-c2c9-11f1-2151-a1386d958103
# ╟─32de4ac0-c2c9-11f1-132c-811be0ca02c5
# ╟─32de4ac0-c2c9-11f1-1109-413fcf0723ee
# ╟─32de4ac0-c2c9-11f1-1028-e39ac09ff156
# ╟─32de4ac0-c2c9-11f1-3acb-2f01b10eb837
# ╟─32de4ac0-c2c9-11f1-0f14-95f84c34976c
# ╟─32de4ac0-c2c9-11f1-3e7b-5d2fca3bf741
# ╟─32de4ac0-c2c9-11f1-1dc7-fbb2043060d6
# ╟─32de4ac0-c2c9-11f1-22a3-89e1202105ea
# ╟─32de4ac0-c2c9-11f1-3910-715c854fa3ba
# ╟─32de4ac0-c2c9-11f1-1eb9-65b38f857339
# ╟─32de4ac0-c2c9-11f1-35a9-8f5fef32f61e
