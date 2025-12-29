module PhantomRevealer_GizmoExtension
# Include the Julia Module
# With the order of level
_module_location = @__DIR__

#Level 1 (Package and information)
include("$_module_location/julia/module_initialization.jl")

const ITP_PRIM  = Ref{Union{Nothing, Interpolations.Extrapolation}}(nothing)
const ITP_METAL = Ref{Union{Nothing, Interpolations.Extrapolation}}(nothing)
const ITP_MMW   = Ref{Union{Nothing, Interpolations.Extrapolation}}(nothing)
#Level 2 (Radiative Cooling table)
include("$_module_location/julia/radiative_cooling.jl")


#Level 3.5 (Extending PhantomRevealerDataStructures)
include("$_module_location/julia/add_quantities.jl")
#Level 4 (Sigal point analysis and File read)
include("$_module_location/julia/read_gizmo.jl")
#Level 5 (Estimation tools)
include("$_module_location/julia/common_tools.jl")

assign_cooling_tables(; table = "$_module_location/table/CloudyData_noUVB.h5")

# Export function, marco, const...
for name in filter(s -> !startswith(string(s), "#"), names(@__MODULE__, all = true))
    if !startswith(String(name), "_") && (name != :eval) && (name != :include)
        @eval export $name
    end
end
end