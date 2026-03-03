module Timescales
using HDF5

# Radiative cooling
include(joinpath(@__DIR__, "radiative_cooling", "cooling_table.jl"))
include(joinpath(@__DIR__, "radiative_cooling", "radiative_cooling.jl"))

# Adiabatic cooling
include(joinpath(@__DIR__, "adiabatic_cooling", "adiabatic_cooling.jl"))

# Thermal sputtering
include(joinpath(@__DIR__, "thermal_sputtering", "thermal_sputtering.jl"))

# Non-thermal sputtering
include(joinpath(@__DIR__, "non-thermal_sputtering", "non-thermal_sputtering.jl"))

# total sputtering
include(joinpath(@__DIR__, "total_sputtering", "total_sputtering.jl"))



# Export function, marco, const...
for name in filter(s -> !startswith(string(s), "#"), names(@__MODULE__, all = true))
    if !startswith(String(name), "_") && (name != :eval) && (name != :include)
        @eval export $name
    end
end
end