# TODO: Revert changes or include local system
inpath = "/Users/ssundar/Library/CloudStorage/OneDrive-NREL/codes/ReEDS_nrelgit/runs/hydupgrade_stress_WECC"

ty = 2025
tp = joinpath(inpath)
WEATHERYEAR = 2008 #2007-2013

plant_cf = DataFrame(CSV.File("/Users/ssundar/Library/CloudStorage/OneDrive-NREL/codes/ReEDS-2.0_Input_Processing/hydro-capacityfactors/plant_cf_inreg.csv"))
rename!(plant_cf,:season => :szn)
rename!(plant_cf,:plant_cf_inreg => :value)
plant_cf.tech = lowercase.(plant_cf.tech)
plant_cf.EIA_PtID = Int.(plant_cf.EIA_PtID)

psys_LD = ReEDS2PRAS.reeds_to_pras(inpath, ty, 8760, WEATHERYEAR,hd_existingplant_cf=plant_cf
)

compare_generator_capacities(psys_LD, tp, ty)
compare_line_capacities(psys_LD, tp, ty)

#Testing PRAS Analytics can be done, but not as part of ReEDS2PRAS
# include("/projects/ntps/llavin/PRAS-Analytics/src/PRAS_Analytics.jl") #will only work if you have this path for this test
# PRAS.savemodel(psys_LD,joinpath(tp,"outputs","v20221201_NTPh0_VSC_DemHi_90by2035EarlyPhaseout__core_"*string(ty)*".pras"))#save the system to a .pras file

# pras_sys_location = joinpath(tp,"outputs","v20221201_NTPh0_VSC_DemHi_90by2035EarlyPhaseout__core_"*string(ty)*".pras")

# PRAS_Analytics.run_pras_analysis(pras_sys_location, "v20221201_NTPh0_VSC_DemHi_90by2035EarlyPhaseout__core", 20232008, 10, plots = true, results_location = tp) 
