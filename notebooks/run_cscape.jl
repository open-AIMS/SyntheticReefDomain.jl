# Launch the reef builder + C~scape (R) export notebook in Pluto (opens a browser tab):
#
#   julia --project=notebooks notebooks/run_cscape.jl
#
# It first copies any changes to the domain builder (reef_builder.jl) into cscape_export.jl
# (see sync_cscape_notebook.jl), so the C~scape notebook always has the latest builder.

using Pluto

try
    include(joinpath(@__DIR__, "sync_cscape_notebook.jl"))
catch err
    @warn "Could not update cscape_export.jl from reef_builder.jl; opening it as it is." exception = err
end

Pluto.run(; notebook=joinpath(@__DIR__, "cscape_export.jl"))
