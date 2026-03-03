"""
    YCnonthermal(x::T) where {T<:AbstractFloat}

Evaluate the non-thermal sputtering coefficient for carbon (C) grains as a function of the relative velocity.

This function implements the polynomial fit adopted from Hu et al. (2019) for the non-thermal sputtering
coefficient of carbon grains. The input `x` is defined as `x = log10(v_rel)`, where `v_rel` is the
gas–grain relative velocity (in the same units as used in Hu et al. (2019) for this fit).

The returned value has units of `μm yr⁻¹ cm³`.

# Parameters
- `x :: T`: Logarithmic relative velocity, defined as `x = log10(v_rel)`.

# Returns
- `T`: Non-thermal sputtering coefficient in units of `μm yr⁻¹ cm³`.
"""
@inline function YCnonthermal(x :: T) where {T<:AbstractFloat}
    a0 = T(-4.87021030e+01)
    a1 = T(5.57114303e+01)
    a2 = T(-2.88275125e+01)
    a3 = T(7.44225761e+00)
    a4 = T(-9.65469006e-01)
    a5 = T(5.00882872e-02)

    p = ((((a5*x + a4)*x + a3)*x + a2)*x + a1)*x + a0
    return exp10(p)
end

"""
    YSinonthermal(x::T) where {T<:AbstractFloat}

Evaluate the non-thermal sputtering coefficient for silicon (Si) grains as a function of the relative velocity.

This function implements the polynomial fit adopted from Hu et al. (2019) for the non-thermal sputtering
coefficient of silicon grains. The input `x` is defined as `x = log10(v_rel)`, where `v_rel` is the
gas–grain relative velocity (in the same units as used in Hu et al. (2019) for this fit).

The returned value has units of `μm yr⁻¹ cm³`.

# Parameters
- `x :: T`: Logarithmic relative velocity, defined as `x = log10(v_rel)`.

# Returns
- `T`: Non-thermal sputtering coefficient in units of `μm yr⁻¹ cm³`.
"""
@inline function YSinonthermal(x :: T) where {T<:AbstractFloat}
    a0 = T(-3.24171434e+01)
    a1 = T(2.86540492e+01)
    a2 = T(-1.15865054e+01)
    a3 = T(2.23338917e+00)
    a4 = T(-2.08259379e-01)
    a5 = T(7.36059545e-03)

    p = ((((a5*x + a4)*x + a3)*x + a2)*x + a1)*x + a0
    return exp10(p)
end

"""
    nonthermal_sputtering_timescale(::Val{:C}, nH::T, vrel::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the non-thermal sputtering timescale (in Myr) for carbon (C) grains.

This function evaluates the relative-velocity-dependent non-thermal sputtering coefficient using
`YCnonthermal(log10(vrel))`, based on the fitting formula adopted from Hu et al. (2019), and computes
the corresponding sputtering timescale for the given hydrogen nuclei number density `nH`, relative
velocity `vrel`, and grain radius `a`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `vrel :: T`: Gas–grain relative velocity in km s⁻¹.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering timescale calculation.

# Returns
- `T`: Non-thermal sputtering timescale in Myr.
"""
@inline function nonthermal_sputtering_timescale(::Val{:C}, nH :: T, vrel :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Ynth = T(1e6) * YCnonthermal(log10(vrel))
    return (0.33 * a) / (nH * Ynth)                 # Myr
end

"""
    specific_nonthermal_sputtering_rate(::Val{:C}, nH::T, vrel::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the specific (fractional) non-thermal sputtering rate for carbon (C) grains.

This function evaluates the relative-velocity-dependent non-thermal sputtering coefficient using
`YCnonthermal(log10(vrel))`, based on the fitting formula adopted from Hu et al. (2019), and converts it
to a grain-radius loss rate per unit grain radius for the given hydrogen nuclei number density `nH`,
relative velocity `vrel`, and grain radius `a`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `vrel :: T`: Gas–grain relative velocity in km s⁻¹.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering rate calculation.

# Returns
- `T`: Specific non-thermal sputtering rate in Myr⁻¹.
"""
@inline function specific_nonthermal_sputtering_rate(::Val{:C}, nH :: T, vrel :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Ynth = T(1e6) * YCnonthermal(log10(vrel))
    return (nH * Ynth) / (0.33 * a)                 # Myr
end

"""
    nonthermal_sputtering_timescale(::Val{:Si}, nH::T, vrel::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the non-thermal sputtering timescale (in Myr) for silicon (Si) grains.

This function evaluates the relative-velocity-dependent non-thermal sputtering coefficient using
`YSinonthermal(log10(vrel))`, based on the fitting formula adopted from Hu et al. (2019), and computes
the corresponding sputtering timescale for the given hydrogen nuclei number density `nH`, relative
velocity `vrel`, and grain radius `a`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `vrel :: T`: Gas–grain relative velocity in km s⁻¹.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering timescale calculation.

# Returns
- `T`: Non-thermal sputtering timescale in Myr.
"""
@inline function nonthermal_sputtering_timescale(::Val{:Si}, nH :: T, vrel :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Ynth = T(1e6) * YSinonthermal(log10(vrel))
    return (0.33 * a) / (nH * Ynth)                 # Myr
end


"""
    specific_nonthermal_sputtering_rate(::Val{:Si}, nH::T, vrel::T; a::T = T(0.035)) where {T<:AbstractFloat}

Return the specific (fractional) non-thermal sputtering rate for silicon (Si) grains.

This function evaluates the relative-velocity-dependent non-thermal sputtering coefficient using
`YSinonthermal(log10(vrel))`, based on the fitting formula adopted from Hu et al. (2019), and converts it
to a grain-radius loss rate per unit grain radius for the given hydrogen nuclei number density `nH`,
relative velocity `vrel`, and grain radius `a`.

# Parameters
- `nH :: T`: Hydrogen nuclei number density in cm⁻³.
- `vrel :: T`: Gas–grain relative velocity in km s⁻¹.

# Keyword Arguments
- `a :: T = T(0.035)`: Grain radius in μm used in the sputtering rate calculation.

# Returns
- `T`: Specific non-thermal sputtering rate in Myr⁻¹.
"""
@inline function specific_nonthermal_sputtering_rate(::Val{:Si}, nH :: T, vrel :: T; a :: T = T(0.035)) where {T <: AbstractFloat}
    Ynth = T(1e6) * YSinonthermal(log10(vrel))
    return (nH * Ynth) / (0.33 * a)                 # Myr
end