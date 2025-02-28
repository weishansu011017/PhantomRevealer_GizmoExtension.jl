function Calculate_SFR(datas::PhantomRevealerDataFrame,time_interval::Float64=0.0001)
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
    SFRave *= 10.22495 # Transfer from (10^10 M⊙/ 0.978 Gyr) to (M⊙/yr)
    return SFRave
end