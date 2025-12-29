


function read_cooling_tables(;TablePath = "./CloudyData_noUVB.h5")
    return h5open(TablePath, "r") do h5
        cooling_prim = read(h5["CoolingRates/Primordial/Cooling"])
        cooling_metal = read(h5["CoolingRates/Metals/Cooling"])
        cooling_MMW = read(h5["CoolingRates/Primordial/MMW"])

        hden_grid = read(HDF5.attributes(h5["CoolingRates/Primordial/Cooling"])["Parameter1"])
        T_grid = read(HDF5.attributes(h5["CoolingRates/Primordial/Cooling"])["Temperature"])

        itp_prim = extrapolate(interpolate((T_grid, hden_grid), cooling_prim, Gridded(Linear())), Flat())
        itp_metal = extrapolate(interpolate((T_grid, hden_grid), cooling_metal, Gridded(Linear())), Flat())
        itp_MMW = extrapolate(interpolate((T_grid, hden_grid), cooling_MMW, Gridded(Linear())), Flat())
        return itp_prim, itp_metal, itp_MMW
    end
end

function assign_cooling_tables(; table = "./CloudyData_noUVB.h5")
    if isnothing(ITP_PRIM[])           
        @info "Reading cooling tables…" table
        itp_prim, itp_metal, itp_MMW = read_cooling_tables(; TablePath = table)
        ITP_PRIM[]  = itp_prim         
        ITP_METAL[] = itp_metal
        ITP_MMW[] = itp_MMW
    end
end