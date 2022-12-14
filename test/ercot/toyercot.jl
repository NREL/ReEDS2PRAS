
#HPC path for a ReEDS case I loaded in. Eventually start to move to test in this folder
# testpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/_ercot_seq")
NEMS_path = joinpath("/projects/ntps/llavin/ReEDS-2.0")
testpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/erc_conv_2_ercot_seq")
# fpath = joinpath(testpath,"ReEDS_Augur","augur_data")
test_year = 2028

ReEDS2PRAS.make_pras_system_from_mapping_info(testpath,test_year,"ERCOT_"*string(test_year),NEMS_path)
# ReEDS2PRAS.make_pras_system_from_mapping_info(testpath,2029,"ERCOT_2029")