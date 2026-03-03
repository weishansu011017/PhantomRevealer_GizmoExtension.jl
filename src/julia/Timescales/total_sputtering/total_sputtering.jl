"""
    sputtering_timescale(::Val{:C}, nH::T, Tg::T, vrel::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the total (thermal + non-thermal) sputtering timescale (in Myr) for carbon (C) grains.

This function evaluates the thermal sputtering coefficient using `YCthermal(log10(Tg))` and the non-thermal
sputtering coefficient using `YCnonthermal(log10(vrel))`, based on the fitting formulae adopted from
Hu et al. (2019). The total sputtering coefficient is computed as `Ytot = Yth + Ynth`, and the timescale
is then evaluated for the given hydrogen nuclei number density `nH`, gas temperature `Tg`, relative velocity
`vrel`, and grain radius `a`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `Tg :: T`: Gas temperature in Kelvin.
- `vrel :: T`: Gas–grain relative velocity in km s⁻¹.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering timescale calculation.

# Returns
- `T`: Total sputtering timescale in Myr.
"""
@inline function sputtering_timescale(::Val{:C}, nH :: T, Tg :: T, vrel :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Yth = T(1e6) * YCthermal(log10(Tg))
    Ynth = T(1e6) * YCnonthermal(log10(vrel))
    Ytot = Yth + Ynth
    return (0.33 * a) / (nH * Ytot)                 # Myr
end

"""
    sputtering_timescale(::Val{:Si}, nH::T, Tg::T, vrel::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the total (thermal + non-thermal) sputtering timescale (in Myr) for silicon (Si) grains.

This function evaluates the thermal sputtering coefficient using `YSithermal(log10(Tg))` and the non-thermal
sputtering coefficient using `YSinonthermal(log10(vrel))`, based on the fitting formulae adopted from
Hu et al. (2019). The total sputtering coefficient is computed as `Ytot = Yth + Ynth`, and the timescale
is then evaluated for the given hydrogen nuclei number density `nH`, gas temperature `Tg`, relative velocity
`vrel`, and grain radius `a`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `Tg :: T`: Gas temperature in Kelvin.
- `vrel :: T`: Gas–grain relative velocity in km s⁻¹.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering timescale calculation.

# Returns
- `T`: Total sputtering timescale in Myr.
"""
@inline function sputtering_timescale(::Val{:Si}, nH :: T, Tg :: T, vrel :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Yth = T(1e6) * YSithermal(log10(Tg))
    Ynth = T(1e6) * YSinonthermal(log10(vrel))
    Ytot = Yth + Ynth
    return (0.33 * a) / (nH * Ytot)                 # Myr
end