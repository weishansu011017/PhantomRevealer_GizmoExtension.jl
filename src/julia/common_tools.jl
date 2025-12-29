function Calculate_SFR(datas::ParticleDataFrame,time_interval::Float64=0.0001)
    # Output would be in (M⊙/yr)
    current_time = datas.params["Time"]
    ntime = Int64(floor(current_time/time_interval))
    SFR_array = zeros(Float64,ntime)
    Star_formation_time = datas[:,"StellarFormationTime"]
    Star_mass = datas[:,"m"]
    for i in eachindex(SFR_array)
        SFR_array[i] = sum(Star_mass[(Star_formation_time .> ((i-1)*time_interval)) .& (Star_formation_time .< (i*time_interval))]) / time_interval
    end
    SFRave = mean(SFR_array)

    utime = datas.params["utime"]
    umass = datas.params["umass"]
    gls2Msunlyr = 1.587077215067404e-26
    code2Msunlyr = (umass/utime) * gls2Msunlyr

    SFRave *= code2Msunlyr # Transfer from (10^10 M⊙/ 0.978 Gyr) to (M⊙/yr)
    return SFRave
end

function Calculate_cooling_timescale()
end