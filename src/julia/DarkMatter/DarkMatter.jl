"""
A toolbox for dark matter estimation
    Wei-Shan Su, 2026
"""

module DarkMatter

# HernquistProfile
include(joinpath(@__DIR__, "HernquistProfile.jl"))

# Export function, marco, const...
for name in filter(s -> !startswith(string(s), "#"), names(@__MODULE__, all = true))
    if !startswith(String(name), "_") && (name != :eval) && (name != :include)
        @eval export $name
    end
end
end