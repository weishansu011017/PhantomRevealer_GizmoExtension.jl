"""
    Adding extra quantites to the given `PhantomRevealerDataFrame`
        By Wei-Shan Su, 2025
"""

"""
    add_SpeciesGasMassDensity!(data::PhantomRevealerDataFrame; fH :: T = 0.7381, fHe :: T = 0.2485) where {T<:AbstractFloat}
Add the mass and density of H, He for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 

# Keyword Arguments
- `fH :: T = 0.7381`: Desired mass fraction of Hydrogen. Defaults to the solar value (X = 0.7381 (Asplund(2009))).
- `fHe :: T = 0.2485`: Desired mass fraction of Helium. Defaults to the solar value (Y = 0.2485 (Asplund(2009))).

"""
function add_SpeciesGasMassDensity!(data::PhantomRevealerDataFrame; fH :: T = 0.7381, fHe :: T = 0.2485) where {T<:AbstractFloat}
    if (data.params["PartType"] == "PartType0")
        m = data[!,"m"]
        ρ = data[!,"rho"]
        Ftype = eltype(m)
        TfH = Ftype(fH)
        TfHe = Ftype(fHe)
        mH    = similar(m)
        mHe   = similar(m)
        rhoH  = similar(m)
        rhoHe = similar(m)
        @threads for i in eachindex(m)
            @inbounds begin
                mi = m[i]
                ρi = ρ[i]
                mH[i]    = TfH    * mi
                mHe[i]   = TfHe   * mi
                rhoH[i]  = TfH    * ρi
                rhoHe[i] = TfHe   * ρi
            end
        end
        data[!,"mH"] = mH
        data[!,"mHe"] = mHe
        data[!,"rhoH"] = rhoH
        data[!,"rhoHe"] = rhoHe

        data.params["fH"] = fH
        data.params["fHe"] = fHe
    end
end

"""
    add_NumberDensityGas!(data::PhantomRevealerDataFrame)

Compute and add the number densities of hydrogen-bearing species (HI, HII, H₂), helium (He), and carbon monoxide (CO) for each gas particle.

**Note**: Abundances χH2, χH+, χCO are interpreted as hydrogen-nucleus-based values (e.g., χH2 = nH2 / nH).

# Parameters
- `data::PhantomRevealerDataFrame`: The SPH data stored in a PhantomRevealerDataFrame.
"""
function add_NumberDensityGas!(data::PhantomRevealerDataFrame)
    if (data.params["PartType"] == "PartType0")
        if !hasproperty(data.dfdata, "rhoH")
            add_SpeciesGasMassDensity!(data)
        end
        rhoH = data[!,"rhoH"]
        rhoHe = data[!,"rhoHe"]
        Ftype = eltype(rhoH)
        umass = data.params["umass"]
        invmp = umass/1.67262192e-24
        inv4mp = invmp * 0.25
        χH2 = data[!,"χH2"]
        χHp = data[!,"χH+"]
        χCO = data[!,"χCO"]

        nH  = zeros(Float64, length(rhoH))              # Eazy to be overflowed, use Float64
        nHe = similar(nH)
        nHp = similar(nH)
        nH2 = similar(nH)
        nHI = similar(nH)
        nCO = similar(nH)
        χHI = similar(nH)

        @threads for i in eachindex(nH)
            @inbounds begin
                rhoHi = rhoH[i]
                rhoHei = rhoHe[i]
                χHpi = χHp[i]
                χH2i = χH2[i]
                χCOi = χCO[i]
                χHIi = 1 - 2*χH2i - χHpi
                nHi = rhoHi * invmp

                nH[i] = nHi
                nHe[i] = rhoHei * inv4mp
                nHp[i] = nHi * χHpi
                nH2[i] = nHi * χH2i
                nHI[i] = nHi * χHIi
                nCO[i] = nHi * χCOi
                χHI[i] = χHIi
            end
        end

        data[!,"nH"] = nH
        data[!,"nHe"] = nHe
        data[!,"nH+"] = nHp
        data[!,"nH2"] = nH2
        data[!,"nHI"] = nHI
        data[!,"nCO"] = nCO
        data[!,"χHI"] = χHI
    end
end

"""
    add_MeanMolecularWeight!(data::PhantomRevealerDataFrame)
Add the mean molecular weight gas particles for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 

"""
function add_MeanMolecularWeight!(data::PhantomRevealerDataFrame)
    if (data.params["PartType"] == "PartType0")
        if !hasproperty(data.dfdata, "nH")
            add_NumberDensityGas!(data)
        end
        Ftype = eltype(data[!,"nH"])
        # Number densities
        nHI = data[!,"nHI"]
        nHp = data[!,"nH+"]
        nH2 = data[!,"nH2"]
        nHe = data[!,"nHe"]
        nCO = data[!,"nCO"]
        ne = nHe

        # Proton per mass
        mHI = one(Ftype)
        mHp = one(Ftype)
        mH2 = Ftype(2.0)
        mHe = Ftype(4.0)
        mCO = Ftype(28.0)

        μ = similar(nHI)

        @threads for i in eachindex(μ)
            @inbounds begin
                nHIi = nHI[i]
                nHpi = nHp[i]
                nH2i = nH2[i]
                nHei = nHe[i]
                nCOi = nCO[i]
                nei  = ne[i]
                mass = mHI * nHIi + mHp * nHpi + mH2 * nH2i + mHe * nHei + mCO * nCOi 
                numprotons = nHIi + nHpi + nH2i + nHei + nCOi + nei
                μ[i] = mass / numprotons
            end
        end
        data[!,"μ"] = μ
    end
end
"""
    add_GasTemperature!(data::PhantomRevealerDataFrame; adiabatic_index::T = 1.6666666667) where {T<:AbstractFloat}
Add the Gas Temperture for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 

# Keyword arguments
- `adiabatic_index :: T = 1.6666666667`: The adiabatic index of system. Default to be non-relativistic monatomic gas(γ = 5/3).
"""
function add_GasTemperature!(data::PhantomRevealerDataFrame; adiabatic_index::T = 1.6666666667) where {T<:AbstractFloat}
    if (data.params["PartType"] == "PartType0")
        if !hasproperty(data.dfdata, "μ")  
            add_MeanMolecularWeight!(data)
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
    add_DGR!(data::PhantomRevealerDataFrame)
Add the total mass of dust and dust-to-gas ratio for each particles. Valid only if `data` with `data.params["PartType"] == "PartType0` (Gas particles) is provided.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data that is stored in `PhantomRevealerDataFrame` 
"""
function add_DGR!(data::PhantomRevealerDataFrame)
    if (data.params["PartType"] == "PartType0")
        data[!, :md] = data[!, :mC] .+ data[!, :mSi]
        data[!, :DGR] = data[!, :md] ./ data[!, :m]
        data[!, :DGR_C] = data[!, :mC] ./ data[!, :m]
        data[!, :DGR_Si] = data[!, :mSi] ./ data[!, :m]
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
    end
end

"""
    add_adiabatic_constant!(data :: PhantomRevealerDataFrame; adiabatic_index :: Float64 = 1.6666666667)

Compute and add the adiabatic (entropy) constant `Aγ = (γ - 1)uρ^(1 - γ)` for each gas particle, where `u` is the internal energy per unit mass and `ρ` is the density. This constant remains invariant for adiabatic processes and characterizes the specific entropy of each fluid element.

# Parameters
- `data :: PhantomRevealerDataFrame`: The SPH data stored in `PhantomRevealerDataFrame`.

# Keyword Arguments
| Name | Default | Description |
|:------|:----------|:-------------|
| `adiabatic_index` | `1.6666666667` | The adiabatic index γ of the gas. Default corresponds to a non-relativistic monatomic ideal gas (γ = 5/3). |
"""
function add_adiabatic_constant!(data :: PhantomRevealerDataFrame; adiabatic_index :: Float64 = 1.6666666667)
    if (data.params["PartType"] == "PartType0")
        cons = (adiabatic_index - 1)
        ργ = data[!, "rho"].^(-cons)
        data[!,"Aγ"] = @. cons * data[!,"InternalEnergy"] * ργ            # In the unit of umass^(1-γ) * udist ^(3γ-1) * uv
    end
end


function add_VelocityDivergence!(data :: PhantomRevealerDataFrame)
    if (data.params["PartType"] == "PartType0")
        
    end
end

function add_RadiativeCoolingTimescale!(data :: PhantomRevealerDataFrame)
    if (data.params["PartType"] == "PartType0")
    end
end