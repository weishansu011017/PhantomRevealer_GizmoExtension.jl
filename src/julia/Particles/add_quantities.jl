"""
    Adding extra quantites to the given `ParticleDataFrame`
        By Wei-Shan Su, 2025
"""

"""
    add_GasSpeciesDensities!(data::ParticleDataFrame; χHe::Real = 0.079)

Add derived gas-species fields for each particle, including H/He mass densities, several species number densities, and the mean molecular weight `μ`.
Valid only if `data` with `data.params["PartType"] == "PartType0"` (Gas particles) is provided.

This function assumes hydrogen-normalised number fractions (`χa = na / nH`) are available in `data`:
- `χH2` = nH2 / nH
- `χH+` = nH+ / nH
- `χCO` = nCO / nH

Neutral atomic hydrogen is derived from hydrogen nuclei conservation:
`χHI = max(0, 1 - 2χH2 - χH+)`, hence `nHI = nH * χHI`.

Electron number density follows the simulation closure: `ne = nH+`.

Helium abundance is defined by number: `χHe = nHe / nH`.
With `mHe ≈ 4 mp`, the mass densities satisfy `rhoH = rho/(1+4χHe)` and `rhoHe = 4χHe*rhoH`.

Columns added/overwritten:
- `"χHI"`, `"μ"`, `"rhoH"`, `"rhoHe"`, `"nH"`, `"nHe"`, `"nH+"`, `"nH2"`, `"nHI"`, `"nCO"`, `"ne"`.

# Parameters
- `data :: ParticleDataFrame`: The SPH data that is stored in `ParticleDataFrame`. Must contain `"rho"`, `"χH2"`, `"χH+"`, `"χCO"` and `data.params["umass"]`.

# Keyword Arguments
- `χHe :: Real = 0.079`: Helium abundance by number, defined as `χHe = nHe / nH`.
"""
function add_GasSpeciesDensities!(data::ParticleDataFrame; χHe :: Real = 0.079)
    # χHe = nHe / nH
    # In code, ne = nH+ = nHχH+
    # nH = q * (ρ/mp) = (1/(1+4χHe)) * (ρ/mp)
    # ζH = nH/n = 1/(1 + χHe + χH+ - χH2)
    if (data.params["PartType"] == "PartType0")
        ρ = data[!,"rho"]
        TχHe = Float64(χHe)
        q = inv(1.0 + 4.0 * TχHe)           # Keep Float64
        

        umass = data.params["umass"]
        invmp = umass/1.67262192e-24        # Keep Float64
        qlmp = q * invmp                    # Keep Float64
        
        χH2 = data[!,"χH2"]
        χHp = data[!,"χH+"]
        χCO = data[!,"χCO"]

        rhoH  = similar(ρ)
        rhoHe = similar(ρ)
        μ     = similar(ρ)

        nH    = similar(ρ, Float64)
        nHe   = similar(nH)
        nHp   = similar(nH)
        nH2   = similar(nH)
        nHI   = similar(nH)
        nCO   = similar(nH)
        χHI   = similar(nH)
        ζH    = similar(nH)
        @threads for i in eachindex(ρ)
            @inbounds begin
                # Collecting known properties
                ρi       = Float64(ρ[i])
                χHpi     = Float64(χHp[i])
                χH2i     = Float64(χH2[i])
                χCOi     = Float64(χCO[i])

                # Neutral atomic hydrogen (clump with max function)
                χHIi     = neutral_atomic_hydrogen(χH2i, χHpi)
                
                # Density of various gas species
                rhoHi    = q * ρi
                rhoHei   = 4.0 * rhoHi * TχHe

                # Mean molecular weight
                μi       = mean_molecular_weight(χH2i, χHpi, χHIi, TχHe)

                # H-nucleus fraction of total particle number density: ζH = nH / n
                ζHi      = hydrogen_nucleus_fraction(χH2i, χHpi, TχHe)

                # Number densities of various element
                nHi      = qlmp * ρi
                nHei     = nHi * TχHe
                nHpi     = nHi * χHpi
                nH2i     = nHi * χH2i
                nHIi     = nHi * χHIi
                nCOi     = nHi * χCOi
                
                # Assign back to 
                χHI[i]   = χHIi
                ζH[i]    = ζHi
                μ[i]     = μi
                rhoH[i]  = rhoHi
                rhoHe[i] = rhoHei
                nH[i]    = nHi
                nHe[i]   = nHei
                nHp[i]   = nHpi
                nH2[i]   = nH2i
                nHI[i]   = nHIi
                nCO[i]   = nCOi
            end
        end
        data[!,"χHI"]   = χHI
        data[!,"ζH"]    = ζH
        data[!,"μ"]     = μ
        data[!,"rhoH"]  = rhoH
        data[!,"rhoHe"] = rhoHe
        data[!,"nH"]    = nH
        data[!,"nHe"]   = nHe
        data[!,"nH+"]   = nHp
        data[!,"nH2"]   = nH2
        data[!,"nHI"]   = nHI
        data[!,"nCO"]   = nCO
        data[!,"ne"]    = nHp

        data.params["χHe"] = χHe
    end
    return nothing
end

"""
    add_GasTemperature!(data::ParticleDataFrame; adiabatic_index::T = 1.6666666667) where {T<:AbstractFloat}
Add the Gas Temperture for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: ParticleDataFrame`: The SPH data that is stored in `ParticleDataFrame` 

# Keyword arguments
- `adiabatic_index :: T = 1.6666666667`: The adiabatic index of system. Default to be non-relativistic monatomic gas(γ = 5/3).
"""
function add_GasTemperature!(data::ParticleDataFrame; adiabatic_index::T = 1.6666666667) where {T<:AbstractFloat}
    if (data.params["PartType"] == "PartType0")
        if !hasproperty(data.dfdata, "μ")  
            add_GasSpeciesDensities!(data)
        end     
        InternalEnergy = data[!,"InternalEnergy"]
        μ = data[!, "μ"]
        Ftype = eltype(μ)
        kB = Ftype(1.380649e-16)                        # erg/K
        gamma = Ftype(adiabatic_index)        
        mp = Ftype(1.67262192595e-24)                   # g
                     
        Q = Ftype(1e10)*(mp*(gamma-one(Ftype)))/(kB)    # g*K/erg = 1.0e10 * K/((km/s)^2)
        Tg = similar(μ)
        @threads for i in eachindex(Tg)
            @inbounds begin
                InternalEnergyi = InternalEnergy[i]
                μi = μ[i]
                Tg[i] = Q * InternalEnergyi * μi
            end
        end
        data[!,"GasTemperature"] = Tg
    end
end

"""
    add_DGR!(data::ParticleDataFrame)
Add the total mass of dust and dust-to-gas ratio for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: ParticleDataFrame`: The SPH data that is stored in `ParticleDataFrame` 
"""
function add_DGR!(data::ParticleDataFrame)
    if (data.params["PartType"] == "PartType0")
        ndustspecies = data.params["NumDustSpecies"]
        if ndustspecies == 1
            data[!, :md] = data[!, :mC] .+ data[!, :mSi]
            data[!, :DGR] = data[!, :md] ./ data[!, :m]
            data[!, :DGR_C] = data[!, :mC] ./ data[!, :m]
            data[!, :DGR_Si] = data[!, :mSi] ./ data[!, :m]
        else
            mg      = data.dfdata.m
            npart   = length(mg)
            md      = zeros(Float64, npart)
            DGR_C   = zeros(Float64, npart)
            DGR_Si  = zeros(Float64, npart)
            @inbounds for n in 1:ndustspecies
                mCsymbol = Symbol("mC_$(n)")
                mSisymbol = Symbol("mSi_$(n)")
                mCn = data[!, mCsymbol]
                mSin = data[!, mSisymbol]
                @inbounds @threads for i in 1:npart
                    @inbounds begin
                        mCni = mCn[i]
                        mSini = mSin[i]
                        mdni = mCni + mSini

                        md[i] += mdni 
                        DGR_C[i] += mCni
                        DGR_Si[i] += mSini
                    end
                end
            end
            DGR_C ./= mg
            DGR_Si ./= mg
            DGR = DGR_C .+ DGR_Si

            data[!, :md] = md
            data[!, :DGR] = DGR
            data[!, :DGR_C] = DGR_C
            data[!, :DGR_Si] = DGR_Si
        end
    end
end


"""
    add_Pressure_adiabetic!(data::ParticleDataFrame; adiabatic_index::Float64 = 1.6666666667)
Add the pressure of gas particles for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: ParticleDataFrame`: The SPH data that is stored in `ParticleDataFrame` 

# Keyword arguments
- `adiabatic_index :: Float64 = 1.6666666667`: The adiabatic index of system. Default to be non-relativistic monatomic gas(γ = 5/3).
"""
function add_Pressure_adiabetic!(data::ParticleDataFrame; adiabatic_index::Float64 = 1.6666666667)
    if (data.params["PartType"] == "PartType0")
        data[!,"P"] = (adiabatic_index - 1) .* data[!,"rho"] .* data[!,"InternalEnergy"] * 1.0004250399845307   # 10^10 M⊙ kpc^-3 (km/s)^2  => 10^10 M⊙ kpc^-1 Gyr^-2
    end
end

"""
    add_Soundspeed_adiabetic!(data::ParticleDataFrame; adiabatic_index::Float64 = 1.6666666667)
Add the sound of speed of gas particles for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: ParticleDataFrame`: The SPH data that is stored in `ParticleDataFrame` 

# Keyword arguments
- `adiabatic_index :: Float64 = 1.6666666667`: The adiabatic index of system. Default to be non-relativistic monatomic gas(γ = 5/3).
"""
function add_Soundspeed_adiabetic!(data::ParticleDataFrame; adiabatic_index::Float64 = 1.6666666667)
    if (data.params["PartType"] == "PartType0")
        cons = adiabatic_index*(adiabatic_index - 1)
        data[!,"cs"] = sqrt.(cons .* data[!,"InternalEnergy"])  # In the unit of km/s
    end
end

"""
    add_adiabatic_constant!(data :: ParticleDataFrame; adiabatic_index :: Float64 = 1.6666666667)

Compute and add the adiabatic (entropy) constant `Aγ = (γ - 1)uρ^(1 - γ)` for each gas particle, where `u` is the internal energy per unit mass and `ρ` is the density. This constant remains invariant for adiabatic processes and characterizes the specific entropy of each fluid element.

# Parameters
- `data :: ParticleDataFrame`: The SPH data stored in `ParticleDataFrame`.

# Keyword Arguments
| Name | Default | Description |
|:------|:----------|:-------------|
| `adiabatic_index` | `1.6666666667` | The adiabatic index γ of the gas. Default corresponds to a non-relativistic monatomic ideal gas (γ = 5/3). |
"""
function add_adiabatic_constant!(data :: ParticleDataFrame; adiabatic_index :: Float64 = 1.6666666667)
    if (data.params["PartType"] == "PartType0")
        cons = (adiabatic_index - 1)
        ργ = data[!, "rho"].^(-cons)
        data[!,"Aγ"] = @. cons * data[!,"InternalEnergy"] * ργ            # In the unit of umass^(1-γ) * udist ^(3γ-1) * uv
    end
end

function add_RadiativeCoolingTimescale!(data :: ParticleDataFrame)
    if (data.params["PartType"] == "PartType0")
    end
end