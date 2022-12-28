NP = joinpath("/projects/ntps/llavin/ReEDS-2.0")
inpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs")

# ty = 2030
# tp = joinpath(inpath,"stscen_capconv_Mid_Case")
# psys = ReEDS2PRAS.make_pras_system_from_mapping_info(tp,ty,"StandScen_MidCase_"*string(ty),NP)
# compare_generator_capacities(psys,tp,ty)

ty = 2030
tp = joinpath(inpath,"stscen_capconv_Low_Demand_Growth")
psys_LD = ReEDS2PRAS.make_pras_system_from_mapping_info(tp,ty,"StandScen_Low_Demand_Growth_"*string(ty),NP)
compare_generator_capacities(psys_LD,tp,ty)

ty = 2040
tp = joinpath(inpath,"stscen_capconv_Mid_Case")
psys_2040 = ReEDS2PRAS.make_pras_system_from_mapping_info(tp,ty,"StandScen_MidCase_"*string(ty),NP)
compare_generator_capacities(psys_2040,tp,ty)
