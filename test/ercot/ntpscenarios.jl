NP = joinpath("/projects/ntps/llavin/ReEDS-2.0")
inpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs")

# ty = 2030
# tp = joinpath(inpath,"stscen_capconv_Mid_Case")
# psys = ReEDS2PRAS.make_pras_system_from_mapping_info(tp,ty,"StandScen_MidCase_"*string(ty),NP)
# compare_generator_capacities(psys,tp,ty)

ty = 2047
tp = joinpath(inpath,"ntptest_LCC_DemHi_CurrentPolEarlyPhaseout__core")
psys_LD = ReEDS2PRAS.make_pras_system_from_mapping_info(tp,ty,"LCC_DemHi_CurrentPolEarlyPhaseout_"*string(ty),NP)
compare_line_capacities(psys_LD,tp,ty)
# compare_generator_capacities(psys_LD,tp,ty)

# ty = 2047
# tp = joinpath(inpath,"ntp20_VSC_DemHi_CurrentPolEarlyPhaseout__core")
# psys = ReEDS2PRAS.make_pras_system_from_mapping_info(tp,ty,"VSC_DemHi_CurrentPolEarlyPhaseout_"*string(ty),NP)
# compare_line_capacities(psys,tp,ty)