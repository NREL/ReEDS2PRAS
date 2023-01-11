NP = joinpath("/projects/ntps/llavin/ReEDS-2.0")
inpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs")

###
#here is a list of the current options
#ntptest_AC_DemHi_CurrentPolEarlyPhaseout__core
#ntptest_Xlim_DemHi_CurrentPolEarlyPhaseout__core

#ntp20_Xlim_DemLo_CurrentPolEarlyPhaseout__core
#ntp20_VSC_DemLo_CurrentPolEarlyPhaseout__core

#ntp20_VSC_DemLo_CurrentPolEarlyPhaseout__TransCost5x
###

ty = 2047
tp = joinpath(inpath,"ntp20_VSC_DemLo_CurrentPolEarlyPhaseout__TransCost5x")
psys_LD = ReEDS2PRAS.make_pras_system_from_mapping_info(tp,ty,"VSC_DemLo_CurrentPolEarlyPhaseout_"*string(ty),NP)
compare_line_capacities(psys_LD,tp,ty)
compare_generator_capacities(psys_LD,tp,ty)

# ty = 2038
# tp = joinpath(inpath,"ntp20_VSC_DemHi_CurrentPolEarlyPhaseout__core")
# psys = ReEDS2PRAS.make_pras_system_from_mapping_info(tp,ty,"VSC_DemHi_CurrentPolEarlyPhaseout_"*string(ty),NP)
# compare_line_capacities(psys,tp,ty)