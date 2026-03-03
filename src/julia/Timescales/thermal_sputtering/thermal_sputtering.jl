"""
    YCthermal(x::T) where {T<:AbstractFloat}

Evaluate the thermal sputtering yield factor for carbon (C) as a function of `x = log10(Tg)`.

This function implements the polynomial fit adopted from Hu et al. (2019) for the temperature-dependent
thermal sputtering yield factor of carbon grains. The input `x` is expected to be `log10(Tg)` where `Tg`
is the gas temperature in Kelvin.

The returned value has units of `μm yr⁻¹ cm³`.

# Parameters
- `x :: T`: Logarithmic gas temperature, defined as `x = log10(Tg)`.

# Returns
- `T`: Thermal sputtering coefficient in units of `μm yr⁻¹ cm³`.
"""
@inline function YCthermal(x :: T) where {T<:AbstractFloat}
    a0 = T(-2.34333937e+02)
    a1 = T(1.38485732e+02)
    a2 = T(-3.39021615e+01)
    a3 = T(4.17705353e+00)
    a4 = T(-2.58281473e-01)
    a5 = T(6.38827523e-03)

    p = ((((a5*x + a4)*x + a3)*x + a2)*x + a1)*x + a0
    return exp10(p)
end

"""
    YSithermal(x::T) where {T<:AbstractFloat}

Evaluate the thermal sputtering yield factor for silicon (Si) as a function of `x = log10(Tg)`.

This function implements the polynomial fit adopted from Hu et al. (2019) for the temperature-dependent
thermal sputtering yield factor of silicon grains. The input `x` is expected to be `log10(Tg)` where `Tg`
is the gas temperature in Kelvin.

The returned value has units of `μm yr⁻¹ cm³`.

# Parameters
- `x :: T`: Logarithmic gas temperature, defined as `x = log10(Tg)`.

# Returns
- `T`: Thermal sputtering coefficient in units of `μm yr⁻¹ cm³`.
"""
@inline function YSithermal(x :: T) where {T<:AbstractFloat}
    a0 = T(-2.34790500e+02)
    a1 = T(1.33208637e+02)
    a2 = T(-3.13027448e+01)
    a3 = T(3.71345730e+00)
    a4 = T(-2.21823668e-01)
    a5 = T(5.31746427e-03)

    p = ((((a5*x + a4)*x + a3)*x + a2)*x + a1)*x + a0
    return exp10(p)
end

"""
    thermal_sputtering_timescale(::Val{:C}, nH::T, Tg::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the thermal sputtering timescale (in Myr) for carbon (C) grains.

This function evaluates the temperature-dependent thermal sputtering yield factor using `YCthermal(log10(Tg))`,
based on the fitting formula adopted from Hu et al. (2019). The timescale is then computed from the given hydrogen
nuclei number density `nH`, gas temperature `Tg`, and grain size parameter `a`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `Tg :: T`: Gas temperature in Kelvin.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering timescale calculation.

# Returns
- `T`: Thermal sputtering timescale in Myr.
"""
@inline function thermal_sputtering_timescale(::Val{:C}, nH :: T, Tg :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Yth = T(1e6) * YCthermal(log10(Tg))
    return (0.33 * a) / (nH * Yth)                 # Myr
end

"""
    specific_thermal_sputtering_rate(::Val{:C}, nH::T, Tg::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the specific (fractional) thermal sputtering rate for carbon (C) grains.

This function evaluates the temperature-dependent thermal sputtering coefficient using
`YCthermal(log10(Tg))`, based on the fitting formula adopted from Hu et al. (2019), and converts it to
a grain-radius loss rate per unit grain radius for the given hydrogen nuclei number density `nH` and
gas temperature `Tg`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `Tg :: T`: Gas temperature in Kelvin.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering rate calculation.

# Returns
- `T`: Specific thermal sputtering rate in Myr⁻¹.
"""
@inline function specific_thermal_sputtering_rate(::Val{:C}, nH :: T, Tg :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Yth = T(1e6) * YCthermal(log10(Tg))
    return (nH * Yth) / (0.33 * a)                 # Myr
end

"""
    thermal_sputtering_timescale(::Val{:Si}, nH::T, Tg::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the thermal sputtering timescale (in Myr) for silicon (Si) grains.

This function evaluates the temperature-dependent thermal sputtering yield factor using `YSithermal(log10(Tg))`,
based on the fitting formula adopted from Hu et al. (2019). The timescale is then computed from the given hydrogen
nuclei number density `nH`, gas temperature `Tg`, and grain size parameter `a`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `Tg :: T`: Gas temperature in Kelvin.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering timescale calculation.

# Returns
- `T`: Thermal sputtering timescale in Myr.
"""
@inline function thermal_sputtering_timescale(::Val{:Si}, nH :: T, Tg :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Yth = T(1e6) * YSithermal(log10(Tg))
    return (0.33 * a) / (nH * Yth)                 # Myr
end

"""
    specific_thermal_sputtering_rate(::Val{:Si}, nH::T, Tg::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the specific (fractional) thermal sputtering rate for silicon (Si) grains.

This function evaluates the temperature-dependent thermal sputtering coefficient using
`YSithermal(log10(Tg))`, based on the fitting formula adopted from Hu et al. (2019), and converts it to
a grain-radius loss rate per unit grain radius for the given hydrogen nuclei number density `nH` and
gas temperature `Tg`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `Tg :: T`: Gas temperature in Kelvin.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering rate calculation.

# Returns
- `T`: Specific thermal sputtering rate in Myr⁻¹.
"""
@inline function specific_thermal_sputtering_rate(::Val{:Si}, nH :: T, Tg :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Yth = T(1e6) * YSithermal(log10(Tg))
    return (nH * Yth) / (0.33 * a)                 # Myr
end

