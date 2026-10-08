# Launch the C~scape (R) export notebook in Pluto (opens a browser tab):
#
#   julia --project=notebooks notebooks/run_cscape.jl

using Pluto

Pluto.run(; notebook=joinpath(@__DIR__, "cscape_export.jl"))
