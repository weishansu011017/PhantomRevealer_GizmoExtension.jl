"""
    (rc::RadiativeCoolingEfficiency)(nH, Tg)

Evaluate tabulated radiative cooling quantities at hydrogen number density `nH`
and gas temperature `Tg`.

**API:** call as `rc(nH, Tg)`.

**Interpolation coordinates:** internally, the tables are bilinearly interpolated
on the axes `(Tg, log10(nH))` (first axis: temperature; second axis: `log10(nH)`),
with flat extrapolation outside the tabulated domain.

The callable returns the interpolated primordial and metal-line cooling
efficiencies together with the tabulated mean molecular weight (MMW) from the
same table.

# Parameters
- `nH` : Hydrogen number density (cm⁻³).
- `Tg` : Gas temperature (K).

# Returns
- `Λ_prim` : Primordial cooling efficiency (erg cm³ s⁻¹) at `(Tg, nH)`.
- `Λ_metal`: Metal-line cooling efficiency (erg cm³ s⁻¹) at `(Tg, nH)`.
- `μ`      : Tabulated mean molecular weight (dimensionless) at `(Tg, nH)`.
"""
@inline function (rc::RadiativeCoolingEfficiency{TF})(nH :: TF, Tg :: TF) where {TF <: AbstractFloat}
    if nH <= zero(TF) || Tg <= zero(TF)
        return zero(TF), zero(TF), TF(Inf) 
    end

    lognH = log10(nH)
    # Read grid and corresponding axis
    TgGrid     = rc.Tg
    lognHGrid  = rc.lognH
    primordial = rc.primordial
    metal      = rc.metal
    MMW        = rc.MMW

    # Flat extrapolation
    Tg    = clamp(Tg,    TgGrid[1],    TgGrid[end])
    lognH = clamp(lognH, lognHGrid[1], lognHGrid[end])

    i = searchsortedlast(TgGrid, Tg)
    j = searchsortedlast(lognHGrid, lognH)

    NT = length(TgGrid)
    NN = length(lognHGrid)
    i = clamp(i, 1, NT    - 1)
    j = clamp(j, 1, NN - 1)

    @inbounds begin
        T1, T2 = TgGrid[i], TgGrid[i+1]
        n1, n2 = lognHGrid[j], lognHGrid[j+1]

        p11    = primordial[i, j]
        p21    = primordial[i+1, j]
        p12    = primordial[i, j+1]
        p22    = primordial[i+1, j+1]

        m11    = metal[i,   j]
        m21    = metal[i+1, j]
        m12    = metal[i,   j+1]
        m22    = metal[i+1, j+1]

        u11 = MMW[i,   j]
        u21 = MMW[i+1, j]
        u12 = MMW[i,   j+1]
        u22 = MMW[i+1, j+1]
    end

    wT = (Tg    - T1) / (T2 - T1)
    wn = (lognH - n1) / (n2 - n1)

    pitp = (1 - wT) * (1 - wn) * p11 + wT * (1 - wn) * p21 + (1 - wT) * wn * p12 + wT * wn * p22
    mitp = (1 - wT) * (1 - wn) * m11 + wT * (1 - wn) * m21 + (1 - wT) * wn * m12 + wT * wn * m22
    μitp = (1 - wT) * (1 - wn) * u11 + wT * (1 - wn) * u21 + (1 - wT) * wn * u12 + wT * wn * u22

    return pitp, mitp, μitp
end

"""
    cooling_efficiency(rc, nH, Tg)

Function-style wrapper for the callable interface of `RadiativeCoolingEfficiency`.

This is a thin wrapper that forwards to `rc(nH, Tg)` and returns the same triple
`(Λ_prim, Λ_metal, μ)`.

# Parameters
- `rc` : Radiative cooling table container (callable as `rc(nH, Tg)`).
- `nH` : Hydrogen number density (cm⁻³).
- `Tg` : Gas temperature (K).

# Returns
- `Λ_prim` : Primordial cooling efficiency (erg cm³ s⁻¹) at `(Tg, nH)`.
- `Λ_metal`: Metal-line cooling efficiency (erg cm³ s⁻¹) at `(Tg, nH)`.
- `μ`      : Tabulated mean molecular weight (dimensionless) at `(Tg, nH)`.
"""
@inline cooling_efficiency(rc::RadiativeCoolingEfficiency, nH :: TF, Tg :: TF) where {TF <: AbstractFloat} = rc(nH, Tg)

"""
    radiative_cooling_timescale(rc::RadiativeCoolingEfficiency,
                                nH::TF,
                                Tg::TF;
                                γ::TF = TF(5//3),
                                Z::TF = TF(0.01295)) where {TF<:AbstractFloat}

Return the radiative cooling timescale (in Myr) evaluated from an interpolated cooling efficiency table.

This function queries `rc(nH, Tg)` to obtain the primordial (`pitp`) and metal (`mitp`) cooling efficiencies
and the mean molecular weight `μitp`, then constructs the total cooling efficiency
`ℰitp = pitp + Z * mitp`. The electron abundance is inferred from `μitp` assuming a fixed helium abundance
`χHe = 0.1`, and is clamped to a physically reasonable range. The returned timescale corresponds to the
characteristic time for radiative cooling at the given temperature and hydrogen nuclei number density.

Valid for `Tg > 0` and `nH ≥ 0`. The result becomes `Inf` if the inferred electron density is zero.

# Parameters
- `rc :: RadiativeCoolingEfficiency`: Callable object that returns `(pitp, mitp, μitp)` given `(nH, Tg)`.
- `nH :: TF`: Hydrogen nuclei number density.
- `Tg :: TF`: Gas temperature.

# Keyword Arguments
- `γ :: TF = TF(5//3)`: Adiabatic index used in the definition of the internal-energy timescale.
- `Z :: TF = TF(0.01295)`: Metallicity scaling applied to the metal cooling term returned by `rc`.

# Returns
- `TF`: Radiative cooling timescale in Myr.
"""
function radiative_cooling_timescale(rc::RadiativeCoolingEfficiency,
                                     nH::TF,
                                     Tg::TF;
                                     γ::TF = TF(5//3),
                                     Z::TF = TF(0.01295)) where {TF <: AbstractFloat}
    pitp, mitp, μitp = rc(nH, Tg)
    ℰitp = pitp + Z * mitp                                                                      # Radiative cooling efficiency (erg cm³ s⁻¹)

    s2Myr :: TF = TF(3.168808781402895e-14)
    kB :: TF = TF(1.380649e-16)
    χHe :: TF = TF(0.1)                                                                         # From the definition X ≈ 0.716 in grackle

    χe = (one(TF) + TF(4)*χHe) / μitp - (one(TF) + χHe)
    χe = clamp(χe, zero(TF), one(TF) + TF(2)*χHe)                                               # [0, 1.2] for χHe=0.1
    ne = nH * χe                                                                             
    ζH = μitp * inv(one(TF) + TF(4)*χHe)                                                

    tradcool = s2Myr * (kB * Tg) * inv((γ - one(TF)) * ζH * ne * ℰitp)                         
    return tradcool
end


"""
    specific_radiative_cooling_rate(rc::RadiativeCoolingEfficiency,
                                    nH::TF,
                                    Tg::TF;
                                    γ::TF = TF(5//3),
                                    Z::TF = TF(0.01295)) where {TF<:AbstractFloat}

Return the specific (fractional) radiative cooling rate (in Myr⁻¹) evaluated from an interpolated cooling efficiency table.

This function uses `rc(nH, Tg)` to obtain primordial and metal cooling efficiencies and the mean molecular weight,
then constructs the total cooling efficiency `ℰitp = pitp + Z * mitp`. The electron abundance is inferred from `μitp`
assuming a fixed helium abundance `χHe = 0.1`, and is clamped to a physically reasonable range. The returned value is
the fractional cooling rate corresponding to the inverse of the radiative cooling timescale returned by
`radiative_cooling_timescale` for the same inputs.

Valid for `Tg > 0` and `nH ≥ 0`. The result is `0` if the inferred electron density is zero.

# Parameters
- `rc :: RadiativeCoolingEfficiency`: Callable object that returns `(pitp, mitp, μitp)` given `(nH, Tg)`.
- `nH :: TF`: Hydrogen nuclei number density.
- `Tg :: TF`: Gas temperature.

# Keyword Arguments
- `γ :: TF = TF(5//3)`: Adiabatic index used in the definition of the internal-energy timescale.
- `Z :: TF = TF(0.01295)`: Metallicity scaling applied to the metal cooling term returned by `rc`.

# Returns
- `TF`: Fractional radiative cooling rate in Myr⁻¹.
"""
function specific_radiative_cooling_rate(rc::RadiativeCoolingEfficiency,
                                         nH::TF,
                                         Tg::TF;
                                         γ::TF = TF(5//3),
                                         Z::TF = TF(0.01295)) where {TF <: AbstractFloat}
    pitp, mitp, μitp = rc(nH, Tg)
    ℰitp = pitp + Z * mitp                                                                      # Radiative cooling efficiency (erg cm³ s⁻¹)

    s2Myr :: TF = TF(3.168808781402895e-14)
    kB :: TF = TF(1.380649e-16)
    χHe :: TF = TF(0.1)                                                                         # From the definition X ≈ 0.716 in grackle

    χe = (one(TF) + TF(4)*χHe) / μitp - (one(TF) + χHe)
    χe = clamp(χe, zero(TF), one(TF) + TF(2)*χHe)                                               # [0, 1.2] for χHe=0.1
    ne = nH * χe                                                                             
    ζH = μitp * inv(one(TF) + TF(4)*χHe)                                                              
    
    sradcool = ((γ - one(TF)) * ζH * ne * ℰitp) * inv(kB * Tg * s2Myr)                         
    return sradcool
end


