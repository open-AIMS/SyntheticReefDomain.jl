# `cscape_export.jl` is the domain builder of `reef_builder.jl` plus a C~scape export
# section, and it writes only the C~scape project. This script copies the builder's cells
# into it, so make builder changes in `reef_builder.jl` and then run
#
#   julia --project=notebooks notebooks/sync_cscape_notebook.jl
#
# (`run_cscape.jl`, and so `start_cscape_export.bat`, runs this before opening the notebook.)
#
# Copied: the builder's cells after its title and before its "## Write" heading (so neither
# the builder's title nor its Write section, which writes the ADRIA domain). Kept:
# `cscape_export.jl`'s own title cell (first) and everything from its "## C~scape (R) export"
# heading onwards. All cells in between are replaced.

using Pluto

const BUILDER = joinpath(@__DIR__, "reef_builder.jl")
const COMBINED = joinpath(@__DIR__, "cscape_export.jl")

heading(cell, pattern) = occursin(Regex("^md\"\"\"\\s*" * pattern), lstrip(cell.code))
is_title(cell) = heading(cell, "#\\s")
is_write_heading(cell) = heading(cell, "##\\s+Write\\b")
is_cscape_heading(cell) = heading(cell, "##\\s+C~scape")
first_line(cell) = first(split(strip(cell.code), '\n'))

"Ids of the builder cells that belong in the combined notebook, in order."
function builder_part(nb)
    cells = [nb.cells_dict[id] for id in nb.cell_order]
    stop = something(findfirst(is_write_heading, cells), length(cells) + 1)
    return [c.cell_id for c in cells[1:(stop - 1)] if !is_title(c)]
end

builder = Pluto.load_notebook_nobackup(BUILDER)
combined = Pluto.load_notebook_nobackup(COMBINED)

order = combined.cell_order
start = findfirst(id -> is_cscape_heading(combined.cells_dict[id]), order)
start === nothing && error("$(basename(COMBINED)) has no \"## C~scape (R) export\" heading cell.")
head = [id for id in order[1:(start - 1)] if is_title(combined.cells_dict[id])]
tail = order[start:end]
builder_ids = builder_part(builder)

# Cells of the combined notebook's builder part that are no longer copied
for id in order[1:(start - 1)]
    id in head || id in builder_ids ||
        println("Removed from $(basename(COMBINED)): ", first_line(combined.cells_dict[id]))
end

# A new Notebook, so its cell topology (which sets the order cells are saved in) is rebuilt
cells = [[combined.cells_dict[id] for id in head]; [builder.cells_dict[id] for id in builder_ids];
         [combined.cells_dict[id] for id in tail]]
synced = Pluto.Notebook(cells, COMBINED)
synced.metadata = combined.metadata
synced.nbpkg_ctx = nothing  # the notebook activates its own environment
# Store cells in the order Pluto itself saves them (package cells first), so opening the
# notebook in Pluto does not reorder the file again
synced.topology = Pluto.updated_topology(synced.topology, synced, synced.cells)
synced._cached_topological_order = Pluto.topological_order(synced.topology)

new = sprint(io -> Pluto.save_notebook(io, synced))
if new == read(COMBINED, String)
    println("$(basename(COMBINED)) is already in sync with $(basename(BUILDER)).")
else
    write(COMBINED, new)
    println("Updated $(basename(COMBINED)): $(length(head)) title cell, ",
        "$(length(builder_ids)) builder cells, $(length(tail)) C~scape cells.")
end
