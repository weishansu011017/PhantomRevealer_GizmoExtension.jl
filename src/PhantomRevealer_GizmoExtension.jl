module PhantomRevealer_GizmoExtension
# Include the Julia Module
# With the order of level
_module_location = @__DIR__

#Level 1 (Package and information)
include("$_module_location/julia/module_initialization.jl")
#Level 3.5 (Extending PhantomRevealerDataStructures)
include("$_module_location/julia/add_quantities.jl")
#Level 4 (Sigal point analysis and File read)
include("$_module_location/julia/read_gizmo.jl")
#Level 5 (Estimation tools)
include("$_module_location/julia/common_tools.jl")

# Export function, marco, const...
for name in filter(s -> !startswith(string(s), "#"), names(@__MODULE__, all = true))
    if !startswith(String(name), "_") && (name != :eval) && (name != :include)
        @eval export $name
    end
end
end