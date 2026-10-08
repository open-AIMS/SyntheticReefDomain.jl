# Launch the reef builder notebook in Pluto (opens a browser tab):
#
#   julia --project=notebooks notebooks/run.jl
#
# A script avoids shell quoting differences (Windows PowerShell 5.1 strips the inner double
# quotes from `julia -e '...("...")'`).

using Pluto

Pluto.run(; notebook=joinpath(@__DIR__, "reef_builder.jl"))
