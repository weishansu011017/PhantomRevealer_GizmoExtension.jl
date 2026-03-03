
"""
    RadiativeCoolingEfficiency{TF, VF, MF}

Container for tabulated radiative cooling efficiencies as functions of
gas temperature and hydrogen number density.

This struct stores precomputed cooling tables (primordial and metal-line)
defined on a regular grid in temperature and log hydrogen number density.
The tables are intended to be queried via bilinear interpolation using
the callable interface.

# Parameters
- `TF` : floating-point type used for all stored values.
- `VF` : vector type for grid axes (`Tg`, `lognH`).
- `MF` : matrix type for cooling tables.

# Fields
- `Tg`         : temperature grid.
- `lognH`      : log hydrogen number density grid.
- `primordial` : primordial cooling efficiency table.
- `metal`      : metal-line cooling efficiency table.
- `MMW`        : tabulated mean molecular weight (diagnostic only).
"""
struct RadiativeCoolingEfficiency{TF <: AbstractFloat, VF <: AbstractVector{TF}, MF <: AbstractMatrix{TF}}
    Tg         :: VF
    lognH      :: VF

    primordial :: MF            # size = (length(Tg), length(lognH))
    metal      :: MF            # size = (length(Tg), length(lognH))
    MMW        :: MF            # size = (length(Tg), length(lognH))
end

"""
    read_cooling_efficiency(; TablePath)

Load radiative cooling efficiency tables from an HDF5 file and construct
a `RadiativeCoolingEfficiency` object.

The function reads temperature and density grids from dataset attributes
and loads the corresponding primordial cooling, metal-line cooling, and
tabulated mean molecular weight tables.

# Keyword Arguments
- `TablePath` : path to the HDF5 cooling table file.

# Returns
- `rc` : a `RadiativeCoolingEfficiency` instance populated with table data.
"""
function read_cooling_efficiency(:: Type{T} = Float64 ;TablePath = joinpath(@__DIR__, "../../../table/CloudyData_noUVB.h5")) where {T <: AbstractFloat}
    return HDF5.h5open(TablePath, "r") do h5
        primordial = T.(read(h5["CoolingRates/Primordial/Cooling"]))
        metal      = T.(read(h5["CoolingRates/Metals/Cooling"]))
        MMW        = T.(read(h5["CoolingRates/Primordial/MMW"]))

        lognH      = T.(read(HDF5.attributes(h5["CoolingRates/Primordial/Cooling"])["Parameter1"]))
        Tg         = T.(read(HDF5.attributes(h5["CoolingRates/Primordial/Cooling"])["Temperature"]))
        return RadiativeCoolingEfficiency(Tg, lognH, primordial, metal, MMW)
    end
end