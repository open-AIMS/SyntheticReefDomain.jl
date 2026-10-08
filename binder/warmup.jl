# Opens both notebooks once during the Binder build so the first launch is faster
# (packages, plots and Pluto itself are compiled ahead of time).
#   julia --project=notebooks binder/warmup.jl
import Pluto

session = Pluto.ServerSession()
session.options.server.disable_writing_notebook_files = true
session.options.server.launch_browser = false
for name in ("reef_builder.jl", "cscape_export.jl")
    path = joinpath(@__DIR__, "..", "notebooks", name)
    notebook = Pluto.SessionActions.open(session, path; run_async=false)
    errored = count(c -> c.errored, notebook.cells)
    println("warm-up $name: $(length(notebook.cells)) cells, $errored with errors")
    Pluto.SessionActions.shutdown(session, notebook; async=false)
end
