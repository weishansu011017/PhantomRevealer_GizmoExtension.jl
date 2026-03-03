module IO

using .Threads 
using DataFrames
using HDF5
using PhantomRevealer

# IO for data structure
## Read & Write GIZMO Binary dumpfiles
include(joinpath(@__DIR__, "gizmoIO", "read_gizmo.jl"))


# Export function, marco, const...
for name in filter(s -> !startswith(string(s), "#"), names(@__MODULE__, all = true))
    if !startswith(String(name), "_") && (name != :eval) && (name != :include)
        @eval export $name
    end
end
end