

NP = joinpath("/projects/ntps/llavin/ReEDS-2.0")
inpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs")

ty = 2023
tp = joinpath(inpath,"ntpsrerun_Xlim_DemHi_90by2035EarlyPhaseout__core")
WEATHERYEAR = 2008 #2007-2013
psys_LD = ReEDS2PRAS.reeds_to_pras(tp,ty,NP,8760,WEATHERYEAR)

compare_generator_capacities(psys_LD,tp,ty)
compare_line_capacities(psys_LD,tp,ty)

#Testing PRAS Analytics can be done, but not as part of ReEDS2PRAS
# include("/projects/ntps/llavin/PRAS-Analytics/src/PRAS_Analytics.jl") #will only work if you have this path for this test
# PRAS.savemodel(psys_LD,joinpath(tp,"outputs","v20221201_NTPh0_VSC_DemHi_90by2035EarlyPhaseout__core_"*string(ty)*".pras"))#save the system to a .pras file

# pras_sys_location = joinpath(tp,"outputs","v20221201_NTPh0_VSC_DemHi_90by2035EarlyPhaseout__core_"*string(ty)*".pras")

# PRAS_Analytics.run_pras_analysis(pras_sys_location, "v20221201_NTPh0_VSC_DemHi_90by2035EarlyPhaseout__core", 20232008, 10, plots = true, results_location = tp) 