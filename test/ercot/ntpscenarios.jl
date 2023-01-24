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
tp = joinpath(inpath,"v20221201_NTPh0_Xlim_DemHi_90by2035EarlyPhaseout__core")
psys_LD = ReEDS2PRAS.reeds_to_pras(tp,ty,NP,8760,2008)
# PRAS.savemodel(psys_LD,joinpath(tp,"outputs","VSC_DemLo_CurrentPolEarlyPhaseout_"*string(ty)*".pras"))
short,flow = run_pras_system(psys_LD,10)#just two for now to save time but eventually more
# @info 
compare_line_capacities(psys_LD,tp,ty)
compare_generator_capacities(psys_LD,tp,ty)

# ty = 2038
# tp = joinpath(inpath,"ntp20_VSC_DemHi_CurrentPolEarlyPhaseout__core")
# psys = ReEDS2PRAS.reeds_to_pras(tp,ty,"VSC_DemHi_CurrentPolEarlyPhaseout_"*string(ty),NP)
# compare_line_capacities(psys,tp,ty)