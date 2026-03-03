module Chemistry

# Estimation of chemical abundance
## Hydrogen species
include(joinpath(@__DIR__, "hydrogen_species.jl"))


## Mean molecular weight
include(joinpath(@__DIR__, "mean_molecular_weight.jl"))


# Export function, marco, const...
for name in filter(s -> !startswith(string(s), "#"), names(@__MODULE__, all = true))
    if !startswith(String(name), "_") && (name != :eval) && (name != :include)
        @eval export $name
    end
end
end