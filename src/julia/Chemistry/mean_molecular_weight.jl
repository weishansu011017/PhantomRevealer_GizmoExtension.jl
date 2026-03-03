"""
    mean_molecular_weight(χH2::T, χHp::T, χHI::T, χHe::T) where {T<:AbstractFloat}

Return the mean molecular weight `μ` (in units of `mH`) for a gas mixture
parameterized by hydrogen-normalized number fractions `χa ≡ na / nH`.

Assumptions:
- Electrons come only from hydrogen ionization, `ne = nH+` (so `χe = χHp`).
- Helium is neutral with abundance `χHe ≡ nHe / nH`.

The total particle number density is

`n = nHI + nH2 + nH+ + ne + nHe`,

so

`n / nH = χHI + χH2 + 2χHp + χHe`.

The mass density can be written as

`ρ = mH nH (1 + 4χHe)`,

therefore

`μ = ρ / (mH n) = (1 + 4χHe) / (χHI + χH2 + 2χHp + χHe)`.

# Parameters
- `χH2::T`: Molecular hydrogen fraction, `nH2 / nH`.
- `χHp::T`: Ionized hydrogen fraction, `nH+ / nH`.
- `χHI::T`: Neutral atomic hydrogen fraction, `nHI / nH`.
- `χHe::T`: Helium nuclei fraction, `nHe / nH`.

# Returns
- `μ::T`: Mean molecular weight in units of `mH`.
"""
@inline function mean_molecular_weight(χH2 :: T, χHp :: T, χHI :: T, χHe :: T) where {T <: AbstractFloat}
    invq = one(T) + T(4) * χHe
    denμ = χHI + χH2 + T(2) * χHp + χHe
    return invq / denμ
end

"""
    mean_molecular_weight(χH2::T, χHp::T, χHe::T) where {T<:AbstractFloat}

Return the mean molecular weight `μ` (in units of `mH`) computed from `χH2`, `χHp`,
and `χHe`, where the neutral atomic fraction `χHI` is inferred from H-nuclei
conservation and then used in `mean_molecular_weight(χH2, χHp, χHI, χHe)`.

This method computes

`χHI = 1 - 2χH2 - χHp`

and applies the same non-negative clamping rule as `neutral_atomic_hydrogen`.
As a numerical safeguard, clamping may break exact H-nuclei conservation.

# Parameters
- `χH2::T`: Molecular hydrogen fraction, `nH2 / nH`.
- `χHp::T`: Ionized hydrogen fraction, `nH+ / nH`.
- `χHe::T`: Helium nuclei fraction, `nHe / nH`.

# Returns
- `μ::T`: Mean molecular weight in units of `mH`.
"""
@inline function mean_molecular_weight(χH2 :: T, χHp :: T, χHe :: T) where {T <: AbstractFloat}
    χHI = neutral_atomic_hydrogen(χH2, χHp)
    return mean_molecular_weight(χH2, χHp, χHI, χHe)
end