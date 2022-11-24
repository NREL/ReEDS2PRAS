
#HPC path for a ReEDS case I loaded in. Eventually start to move to test in this folder
testpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/_ercot_seq")
test_year = 2012

ReEDS2PRAS.make_pras_system_from_mapping_info(testpath,test_year,"ERCOT_"*string(test_year))
ReEDS2PRAS.make_pras_system_from_mapping_info(testpath,2020,"ERCOT_2020")