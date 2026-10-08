# Generate a synthetic domain, load it in ADRIA and run counterfactual scenarios for every
# climate level it contains.
#
#   julia --project=sandbox -t 4 sandbox/01_generate_and_run.jl [config.toml]
#
# The domain is written to domains/<name> (replacing an earlier version of the same domain).

using SyntheticReefDomain
using CairoMakie  # adds the plots to the domain's README.md
using ADRIA
using Printf, Statistics

const REPO = dirname(@__DIR__)
config_path = get(ARGS, 1, joinpath(REPO, "configs", "default.toml"))

synthetic = build_domain(config_path)
display(synthetic)
println()

path = write_domain(synthetic, joinpath(REPO, "domains", synthetic.config.name); force=true)
labels = string.(climate_labels(path))
println("Written to $path (climate labels: $(join(labels, ", ")))")

ADRIA.setup()
ENV["ADRIA_OUTPUT_DIR"] = mktempdir()

dom = ADRIA.load_domain(path, first(labels))
scens = ADRIA.sample_cf(dom, 16)
rs = ADRIA.run_scenarios(dom, scens, labels)

# Mean relative cover (fraction of habitable area) per decade, averaged over reefs, for
# each climate level
rc = Array(ADRIA.metrics.relative_cover(rs))  # timesteps x locations x scenarios
scen_labels = string.(Int.(rs.inputs.RCP))  # stored as numbers
yrs = collect(dom.env_layer_md.timeframe)
println("\nMean relative cover by decade")
@printf("%-8s", "years")
foreach(l -> @printf("%10s", "label $l"), labels)
println()
for start in yrs[1]:10:yrs[end]
    rows = findall(y -> start <= y < start + 10, yrs)
    @printf("%-8s", "$(start)s")
    for l in labels
        cols = findall(==(l), scen_labels)
        @printf("%10.3f", mean(rc[rows, :, cols]))
    end
    println()
end
