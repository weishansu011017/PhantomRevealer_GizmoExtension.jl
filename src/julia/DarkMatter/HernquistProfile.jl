"""
    HernquistProfile{T <: AbstractFloat}

Immutable container describing a spherical dark matter halo following the
Hernquist density profile

    ρ(r) = ρ_s / [(r/r_s) (1 + r/r_s)^3].

The total mass is finite and equals `Mtot`. The virial radius `Rvir` is
defined via the overdensity condition

    ρ̄(<Rvir) = Δc ρ_crit,

where

    ρ_crit = 3 H0^2 / (8πG).

All fields are stored using the same floating-point type `T`, ensuring
type stability and numerical consistency.

Units must be self-consistent. With

    G = 4.30091727003628e-6  (kpc km² s⁻² M⊙⁻¹),

the expected unit system is:

- Length: kpc
- Mass: M⊙
- Velocity: km s⁻¹
- H0: km s⁻¹ kpc⁻¹

# Parameters
- `rs::T`    : Hernquist scale radius (kpc)
- `Mtot::T`  : Total halo mass (M⊙)
- `ρs::T`    : Characteristic density, defined as Mtot / (2π rs³)
- `Rvir::T`  : Virial radius (kpc)

# Returns
- `HernquistProfile{T}` : Fully specified halo model with precomputed `ρs`
  and `Rvir`
"""
struct HernquistProfile{T <: AbstractFloat}
    rs      :: T
    Mtot    :: T
    ρs      :: T
    Rvir    :: T
end

"""
    HernquistProfile(rs::T, Mtot::T; Δc::T = T(200), H0::T = T(0.07)) where {T <: AbstractFloat}

Construct a `HernquistProfile` from the scale radius and total mass,
automatically computing the characteristic density `ρs` and the virial
radius `Rvir`.

The characteristic density is defined as

    ρs = Mtot / (2π rs³),

and the virial radius is obtained from the overdensity condition

    ρ̄(<Rvir) = Δc ρ_crit,

with

    ρ_crit = 3 H0² / (8πG).

Internally, `Rvir` is computed by solving

    R (R + rs)² = 2 G Mtot / (Δc H0²)

using the analytic cubic solution.

All quantities must be supplied in a self-consistent unit system. With

    G = 4.30091727003628e-6  (kpc km² s⁻² M⊙⁻¹),

the expected units are:

- Length: kpc
- Mass: M⊙
- H0: km s⁻¹ kpc⁻¹

For reference:

    70 km s⁻¹ Mpc⁻¹ = 0.07 km s⁻¹ kpc⁻¹

# Parameters
- `rs::T`      : Hernquist scale radius (kpc)
- `Mtot::T`    : Total halo mass (M⊙)

# Keyword Arguments
- `Δc::T`      : Overdensity parameter (default 200)
- `H0::T`      : Hubble parameter in km s⁻¹ kpc⁻¹ (default 0.07)

# Returns
- `HernquistProfile{T}` : Fully initialized halo model with precomputed
  `ρs` and `Rvir`
"""
@inline function HernquistProfile(rs :: T, Mtot :: T; Δc :: T = T(200.0), H0 :: T = T(0.07)) where {T <: AbstractFloat}
    ρs  = Mtot * inv(2π * rs * rs * rs)
    Rvir = Rvirial(rs, Mtot, Δc, H0)
    return HernquistProfile(rs, Mtot, ρs, Rvir)
end

"""
    Rvirial(rs::T, Mtot::T, Δc::T, H0::T) where {T <: AbstractFloat}

Compute the virial radius `Rvir` of a Hernquist dark matter halo by solving

    R (R + r_s)^2 = C,

where

    C = 2 G M_tot / (Δc H0^2).

This corresponds to the overdensity definition

    ρ̄(<Rvir) = Δc ρ_crit,

with

    ρ_crit = 3 H0^2 / (8πG).

The cubic equation is solved analytically using Cardano’s formula, yielding
the unique positive real root.

All quantities must be given in a self-consistent unit system. With the
internal value

    G = 4.30091727003628e-6  (kpc km² s⁻² M⊙⁻¹),

the expected units are:

- length: kpc
- mass: M⊙
- H0: km s⁻¹ kpc⁻¹

# Parameters
- `rs::T`     : Hernquist scale radius (kpc)
- `Mtot::T`   : Total halo mass (M⊙)
- `Δc::T`     : Overdensity parameter (typically 200)
- `H0::T`     : Hubble parameter in km s⁻¹ kpc⁻¹

# Returns
- `Rvir::T` : Virial radius in the same length unit as `rs`
"""
@inline function Rvirial(rs :: T, Mtot :: T, Δc :: T, H0 :: T) where {T <: AbstractFloat}
    two  = one(T) + one(T)
    three = two + one(T) 

    G = T(4.30091727003628e-6)          # G = 4.30091727003628e-6 km² s⁻² kpc M⊙⁻¹
    C = two * G * Mtot * inv(Δc * H0 * H0)

    rs2 = rs * rs
    rs3 = rs * rs * rs
    p = -rs2 / three
    q = - (C + two * rs3 * inv(T(27)))

    nql2 = -q / two
    pl3 = p / three
    Δ = nql2 * nql2 + pl3 * pl3 * pl3
    sqrtΔ = sqrt(Δ)

    return cbrt(nql2 + sqrtΔ) + cbrt(nql2 - sqrtΔ) - two * rs / three
end

@inline (HP :: HernquistProfile{T})(r :: Real) where {T <: AbstractFloat} = density(HP, r)

@inline function density(HP :: HernquistProfile{T}, r :: Real) where {T <: AbstractFloat}
    rs = HP.rs
    ρs = HP.ρs

    rT = T(r)
    r̃ = rT * inv(rs)
    r̃p1 = r̃ + one(T)
    invρ̃ = inv(r̃ * r̃p1 * r̃p1 * r̃p1)
    return ρs * invρ̃
end

@inline function enclosure_mass(HP :: HernquistProfile{T}, r :: Real) where {T <: AbstractFloat}
    rs = HP.rs
    Mtot = HP.Mtot

    rT = T(r)
    rsum = rT + rs
    return Mtot * rT * rT * inv(rsum * rsum) 
end

@inline function gravitational_potential(HP :: HernquistProfile{T}, r :: Real) where {T <: AbstractFloat}
    rs = HP.rs
    Mtot = HP.Mtot
    G = T(4.30091727003628e-6)          # G = 4.30091727003628e-6 km² s⁻² kpc M⊙⁻¹

    rT = T(r)
    rsum = rT + rs
    invrsum = inv(rsum)

    GM = G * Mtot
    return - GM * invrsum
end

@inline function gravitational_acceleration(HP :: HernquistProfile{T}, r :: Real) where {T <: AbstractFloat}
    rs = HP.rs
    Mtot = HP.Mtot
    G = T(4.30091727003628e-6)          # G = 4.30091727003628e-6 km² s⁻² kpc M⊙⁻¹

    rT = T(r)
    rsum = rT + rs
    rsum2 = rsum * rsum
    invrsum2 = inv(rsum2)

    GM = G * Mtot
    return - GM * invrsum2
end

@inline function escape_velocity(HP :: HernquistProfile{T}, r :: Real) where {T <: AbstractFloat}
    rs = HP.rs
    Mtot = HP.Mtot
    G = T(4.30091727003628e-6)          # G = 4.30091727003628e-6 km² s⁻² kpc M⊙⁻¹

    rT = T(r)
    rsum = rT + rs
    invrsum = inv(rsum)

    GM = G * Mtot
    return sqrt(T(2) * GM * invrsum)
end