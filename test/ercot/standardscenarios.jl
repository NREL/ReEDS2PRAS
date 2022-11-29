testpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs")

test_year = 2020
ReEDS2PRAS.make_pras_system_from_mapping_info(joinpath(testpath,"_Mid_Case"),test_year,"StandScen_MidCase_"*string(test_year))
ReEDS2PRAS.make_pras_system_from_mapping_info(joinpath(testpath,"_High_Electrification"),test_year,"StandScen_Low_Demand_Growth_"*string(test_year))

test_year = 2030
ReEDS2PRAS.make_pras_system_from_mapping_info(joinpath(testpath,"_Mid_Case"),test_year,"StandScen_MidCase_"*string(test_year))
# ReEDS2PRAS.make_pras_system_from_mapping_info(testpath,2020,"ERCOT_2020")