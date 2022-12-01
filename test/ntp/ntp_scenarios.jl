testpath = joinpath("/projects/ntps/sdhulipa/ReEDS2PRAS-Setup/ReEDS_Output")

test_year = 2030
ReEDS2PRAS.make_pras_system_from_mapping_info(joinpath(testpath,"v20220627_NTPg0_AC_DemHi_100by2035__core"),test_year,"AC_"*string(test_year))
# ReEDS2PRAS.make_pras_system_from_mapping_info(joinpath(testpath,"_High_Electrification"),test_year,"StandScen_Low_Demand_Growth_"*string(test_year))