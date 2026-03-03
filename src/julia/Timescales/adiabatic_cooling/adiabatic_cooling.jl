"""
    adiabatic_cooling_timescale(divv::T, γ::T = T(5//3)) where {T<:AbstractFloat}

Return the adiabatic cooling timescale for an expanding flow, computed from the local velocity divergence.

This is defined only for expansion (`divv > 0`). For `divv ≤ 0` (no expansion or compression),
the function returns `Inf` to indicate that no adiabatic cooling timescale applies.

# Parameters
- `divv :: T`: Velocity divergence (`∇⋅v`). Should have units of s⁻¹.

# Keyword Arguments
- `γ :: T = T(5//3)`: Adiabatic index.

# Returns
- `T`: Adiabatic cooling timescale with units of Myr. Returns `Inf` when `divv ≤ 0`.
"""
@inline function adiabatic_cooling_timescale(divv :: T, γ :: T = T(5//3)) where {T <: AbstractFloat}
    s2Myr :: T = T(3.168808781402895e-14)
    return divv > 0 ? s2Myr * inv((γ - one(T)) * divv) : T(Inf)
end

"""
    specific_adiabatic_cooling_rate(divv::T, γ::T = T(5//3)) where {T<:AbstractFloat}

Return the specific (fractional) adiabatic cooling rate for an expanding flow, computed from the local velocity divergence.

This returns a positive cooling rate only for expansion (`divv > 0`). For `divv ≤ 0` (no expansion or compression),
the function returns `0` since adiabatic cooling does not apply (compression corresponds to adiabatic heating).

# Parameters
- `divv :: T`: Velocity divergence (`∇⋅v`). Should have units of s⁻¹.

# Keyword Arguments
- `γ :: T = T(5//3)`: Adiabatic index.

# Returns
- `T`: Fractional adiabatic cooling rate with units of Myr⁻¹. Returns `0` when `divv ≤ 0`.
"""
@inline function specific_adiabatic_cooling_rate(divv :: T, γ :: T = T(5//3)) where {T <: AbstractFloat}
    s2Myr :: T = T(3.168808781402895e-14)
    return divv > 0 ? (γ - one(T)) * divv / s2Myr : zero(T)
end 