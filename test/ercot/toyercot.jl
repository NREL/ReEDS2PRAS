
#HPC path for a ReEDS case I loaded in. Eventually start to move to test in this folder
# testpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/_ercot_seq")
NEMS_path = joinpath("/projects/ntps/llavin/ReEDS-2.0")
testpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/erc_conv_2_ercot_seq")
# fpath = joinpath(testpath,"ReEDS_Augur","augur_data")
test_year = 2028

psys = ReEDS2PRAS.reeds_to_pras(testpath,test_year,NEMS_path,8760,2012)
# ReEDS2PRAS.reeds_to_pras(testpath,2029,"ERCOT_2029")

compare_generator_capacities(psys,testpath,test_year)
compare_line_capacities(psys,testpath,test_year)