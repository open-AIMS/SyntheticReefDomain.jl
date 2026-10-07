# Check that a generated domain translates to C~scape inputs via CscapeInterface.
# C~scape reads each reef's `ub_med` (median wave orbital velocity) from the geopackage, and
# reefs sharing an identical (ub_med, depth_med) pair collapse onto one demography point.
#
# Runs in CscapeInterface's own environment (it is not a dependency of this package):
#
#   julia --project=../CscapeInterface.jl sandbox/02_cscape_translate.jl [domain folder]
#
# Generate the domain first, e.g. with sandbox/01_generate_and_run.jl.

using ADRIA, CscapeInterface

const REPO = dirname(@__DIR__)
path = get(ARGS, 1, joinpath(REPO, "domains", "SyntheticReefs"))
isdir(path) || error("No domain at $path; run sandbox/01_generate_and_run.jl first.")

# First climate label available in the domain (DHWs/dhwRCP<label>.nc)
label = first(sort([
    String(m.captures[1]) for m in (match(r"^dhwRCP(\d+)\.nc$", f) for f in readdir(joinpath(path, "DHWs")))
    if m !== nothing
]; by=x -> parse(Int, x)))

dom = ADRIA.load_domain(path, label)
inputs = CscapeInterface.domain_to_cscape_inputs(dom; dhw_idx=1, cyclone_idx=1)  # wave_idx=0
spatial = inputs.spatial_override
show(stdout, MIME("text/plain"), spatial; allrows=true)
println()

pairs = collect(zip(spatial.ub_med, spatial.depth_med))
println("ub_med taken from the geopackage: ", spatial.ub_med ≈ Float64.(dom.loc_data.ub_med))
println("distinct (ub_med, depth_med) pairs: $(length(unique(pairs))) of $(length(pairs))")
println("translated inputs: ", join(keys(inputs), ", "))
