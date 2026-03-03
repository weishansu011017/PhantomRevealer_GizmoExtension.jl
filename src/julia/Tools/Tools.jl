module Tools

using .Threads
using PhantomRevealer

# Adding quantities function
include(joinpath(@__DIR__, "common_tools.jl"))


# Export function, marco, const...
for name in filter(s -> !startswith(string(s), "#"), names(@__MODULE__, all = true))
    if !startswith(String(name), "_") && (name != :eval) && (name != :include)
        @eval export $name
    end
end
end