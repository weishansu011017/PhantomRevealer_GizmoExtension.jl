module PhantomRevealer_GizmoExtension
# Include the Julia Module
using Pkg
using Reexport

# Tools
include(joinpath(@__DIR__, "julia", "Tools", "Tools.jl"))
@reexport using .Tools

# Chemistry-related tools
include(joinpath(@__DIR__, "julia", "Chemistry", "Chemistry.jl"))
@reexport using .Chemistry

# Toolbox for dark matter
include(joinpath(@__DIR__, "julia", "DarkMatter", "DarkMatter.jl"))
@reexport using .DarkMatter

# Read GIZMO (Chia-Yu Hu Specilized)
include(joinpath(@__DIR__, "julia", "IO", "IO.jl"))
@reexport using .IO

# Add extra quantities
include(joinpath(@__DIR__, "julia", "Particles", "Particles.jl"))
@reexport using .Particles

# Timescale estimation
include(joinpath(@__DIR__, "julia", "Timescales", "Timescales.jl"))
@reexport using .Timescales

end