"""
    Adding extra quantites to the given `PhantomRevealerDataFrame`
        By Wei-Shan Su, 2025
"""

"""
    add_GasTemperature!(data::PhantomRevealerDataFrame)
Add the Gas Temperture for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 
"""
function add_GasTemperature!(data::PhantomRevealerDataFrame)
    if (data.params["PartType"] == "PartType0")
        kB = 1.380649e-16                        # erg/K
        gamma = 5/3           
        mp = 1.67262192595e-24                   # g
        μ = 1.0                     
        Q = 1e10*(μ*mp*(gamma-1))/(kB)           # g*K/erg = 1.0e10 * K/((km/s)^2)
        data[!,"GasTemperature"] = Q .* data[!,"InternalEnergy"]
    else
        nothing
    end
end

"""
    add_DGR!(data::PhantomRevealerDataFrame)
Add the dust-to-gas ratio for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 
"""
function add_DGR!(data::PhantomRevealerDataFrame)
    if (data.params["PartType"] == "PartType0")
        data[!,"DGR"] = (data[!,"CarbonDustMass"] .+ data[!,"SilicateDustMass"])./data[!,"m"]
        data[!,"DGR_C"] = data[!,"CarbonDustMass"]./data[!,"m"]
        data[!,"DGR_Si"] = data[!,"SilicateDustMass"]./data[!,"m"]
    else
        nothing
    end
end


"""
    add_NumberDensityGas!(data::PhantomRevealerDataFrame)
Add the number of density of monatomic gas particles for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 
"""
function add_NumberDensityGas!(data::PhantomRevealerDataFrame)
    if (data.params["PartType"] == "PartType0")
        umass = data.params["umass"]
        mp = 1.6e-24/umass
        data[!,"n"] = (data[!,"rho"] ./mp)
    else
        nothing
    end
end

"""
    add_Pressure_adiabetic!(data::PhantomRevealerDataFrame; adiabatic_index::Float64 = 1.6666666667)
Add the pressure of gas particles for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 

# Keyword arguments
- `adiabatic_index :: Float64 = 1.6666666667`: The adiabatic index of system. Default to be non-relativistic monatomic gas(γ = 5/3).
"""
function add_Pressure_adiabetic!(data::PhantomRevealerDataFrame; adiabatic_index::Float64 = 1.6666666667)
    if (data.params["PartType"] == "PartType0")
        data[!,"P"] = (adiabatic_index - 1) .* data[!,"rho"] .* data[!,"InternalEnergy"] * 1.0004250399845307   # 10^10 M⊙ kpc^-3 (km/s)^2  => 10^10 M⊙ kpc^-1 Gyr^-2
    else
        nothing
    end
end

"""
    add_Soundspeed_adiabetic!(data::PhantomRevealerDataFrame; adiabatic_index::Float64 = 1.6666666667)
Add the sound of speed of gas particles for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 

# Keyword arguments
- `adiabatic_index :: Float64 = 1.6666666667`: The adiabatic index of system. Default to be non-relativistic monatomic gas(γ = 5/3).
"""
function add_Soundspeed_adiabetic!(data::PhantomRevealerDataFrame; adiabatic_index::Float64 = 1.6666666667)
    if (data.params["PartType"] == "PartType0")
        cons = adiabatic_index*(adiabatic_index - 1)
        data[!,"cs"] = sqrt.(cons .* data[!,"InternalEnergy"])  # In the unit of km/s
    else
        nothing
    end
end