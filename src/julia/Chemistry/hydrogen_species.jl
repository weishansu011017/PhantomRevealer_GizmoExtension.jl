"""
    neutral_atomic_hydrogen_raw(χH2::T, χHp::T) where {T<:AbstractFloat}

Return the hydrogen-normalized number fraction of neutral atomic hydrogen, `χHI`,
computed from the H-nuclei conservation constraint

`χHI + χHp + 2χH2 = 1`,

without any clamping.

# Parameters
- `χH2::T`: Hydrogen-normalized number fraction of molecular hydrogen, `nH2 / nH`.
- `χHp::T`: Hydrogen-normalized number fraction of ionized hydrogen, `nH+ / nH`.

# Returns
- `χHI_raw::T`: Unclamped neutral atomic hydrogen fraction.
"""
@inline neutral_atomic_hydrogen_raw(χH2 :: T, χHp :: T) where {T <: AbstractFloat} = one(T) - T(2) * χH2 - χHp


"""
    ionized_hydrogen_raw(χH2::T, χHI::T) where {T<:AbstractFloat}

Return the hydrogen-normalized number fraction of ionized hydrogen, `χHp`,
computed from the H-nuclei conservation constraint

`χHI + χHp + 2χH2 = 1`,

without any clamping.

# Parameters
- `χH2::T`: Hydrogen-normalized number fraction of molecular hydrogen, `nH2 / nH`.
- `χHI::T`: Hydrogen-normalized number fraction of neutral atomic hydrogen, `nHI / nH`.

# Returns
- `χHp_raw::T`: Unclamped ionized hydrogen fraction.
"""
@inline ionized_hydrogen_raw(χH2 :: T, χHI :: T) where {T <: AbstractFloat} = one(T) - T(2) * χH2 - χHI

"""
    hydrogen_molecular_raw(χHp::T, χHI::T) where {T<:AbstractFloat}

Return the hydrogen-normalized number fraction of molecular hydrogen, `χH2`,
computed from the H-nuclei conservation constraint

`χHI + χHp + 2χH2 = 1`,

without any clamping.

# Parameters
- `χHp::T`: Hydrogen-normalized number fraction of ionized hydrogen, `nH+ / nH`.
- `χHI::T`: Hydrogen-normalized number fraction of neutral atomic hydrogen, `nHI / nH`.

# Returns
- `χH2_raw::T`: Unclamped molecular hydrogen fraction.
"""
@inline hydrogen_molecular_raw(χHp :: T, χHI :: T) where {T <: AbstractFloat} = (one(T) - χHp - χHI) / T(2)


"""
    neutral_atomic_hydrogen(χH2::T, χHp::T) where {T<:AbstractFloat}

Return the hydrogen-normalized number fraction of neutral atomic hydrogen, `χHI`,
computed from the H-nuclei conservation constraint

`χHI + χHp + 2χH2 = 1`,

and clamped to be non-negative.

This is a numerical safeguard; clamping may break the exact conservation relation.

# Parameters
- `χH2::T`: Hydrogen-normalized number fraction of molecular hydrogen, `nH2 / nH`.
- `χHp::T`: Hydrogen-normalized number fraction of ionized hydrogen, `nH+ / nH`.

# Returns
- `χHI::T`: Non-negative neutral atomic hydrogen fraction.
"""
@inline neutral_atomic_hydrogen(χH2 :: T, χHp :: T) where {T <: AbstractFloat} = _clamp_nonneg(neutral_atomic_hydrogen_raw(χH2, χHp))

"""
    ionized_hydrogen(χH2::T, χHI::T) where {T<:AbstractFloat}

Return the hydrogen-normalized number fraction of ionized hydrogen, `χHp`,
computed from the H-nuclei conservation constraint

`χHI + χHp + 2χH2 = 1`,

and clamped to be non-negative.

This is a numerical safeguard; clamping may break the exact conservation relation.

# Parameters
- `χH2::T`: Hydrogen-normalized number fraction of molecular hydrogen, `nH2 / nH`.
- `χHI::T`: Hydrogen-normalized number fraction of neutral atomic hydrogen, `nHI / nH`.

# Returns
- `χHp::T`: Non-negative ionized hydrogen fraction.
"""
@inline ionized_hydrogen(χH2 :: T, χHI :: T) where {T <: AbstractFloat} = _clamp_nonneg(ionized_hydrogen_raw(χH2, χHI))


"""
    hydrogen_molecular(χHp::T, χHI::T) where {T<:AbstractFloat}

Return the hydrogen-normalized number fraction of molecular hydrogen, `χH2`,
computed from the H-nuclei conservation constraint

`χHI + χHp + 2χH2 = 1`,

and clamped to be non-negative.

This is a numerical safeguard; clamping may break the exact conservation relation.

# Parameters
- `χHp::T`: Hydrogen-normalized number fraction of ionized hydrogen, `nH+ / nH`.
- `χHI::T`: Hydrogen-normalized number fraction of neutral atomic hydrogen, `nHI / nH`.

# Returns
- `χH2::T`: Non-negative molecular hydrogen fraction.
"""
@inline hydrogen_molecular(χHp :: T, χHI :: T) where {T <: AbstractFloat} = _clamp_nonneg(hydrogen_molecular_raw(χHp, χHI))

"""
    hydrogen_nucleus_fraction(χH2::T, χHp::T, χHe::T) where {T<:AbstractFloat}

Return the hydrogen-nucleus fraction of the total particle number density,

`ζH ≡ nH / n`,

under the assumptions used in this project:
- number fractions are hydrogen-normalized, `χa ≡ na / nH`,
- electrons come only from hydrogen ionization, `ne = nH+` (so `χe = χHp`),
- helium is neutral with fixed abundance `χHe ≡ nHe / nH`.

With `n = nHI + nH2 + nH+ + ne + nHe` and H-nuclei conservation
`χHI + χHp + 2χH2 = 1`, one obtains

`ζH = 1 / (1 + χHe + χHp - χH2)`.

# Parameters
- `χH2::T`: Hydrogen-normalized number fraction of molecular hydrogen, `nH2 / nH`.
- `χHp::T`: Hydrogen-normalized number fraction of ionized hydrogen, `nH+ / nH`.
- `χHe::T`: Hydrogen-normalized number fraction of helium nuclei, `nHe / nH`.

# Returns
- `ζH::T`: Hydrogen-nucleus fraction of the total particle number density, `nH / n`.
"""
@inline hydrogen_nucleus_fraction(χH2 :: T, χHp :: T, χHe :: T) where {T <: AbstractFloat} = inv(one(T) + χHe + χHp - χH2)

# Toolbox
@inline _clamp_nonneg(x :: T) where {T<:AbstractFloat} = max(zero(T), x)