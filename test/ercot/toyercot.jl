
#HPC path for a ReEDS case I loaded in. Eventually start to move to test in this folder
testpath = joinpath("/projects/ntps/llavin/ReEDS-2.0/runs/_ercot_seq")
test_year = 2012

ReEDS2PRAS.make_pras_system_from_mapping_info(testpath,test_year)
ReEDS2PRAS.make_pras_system_from_mapping_info(testpath,2020)

# load_info = h5read(joinpath(testpath,"load.h5"), "data") 
#have to figure how to go from this to relevant pras object

#we can get the regions from this
#assign them loads


#then we want the generation capacity in a given year

#for vg, we will need profiles

#for conventional, we will need failure and recovery probabilities

#for 

# example of how Gord tests PLEXOS2PRAS
# @testset "Toy Model" begin

#     testpath = dirname(@__FILE__) * "/"

#     @testset "Pre-PLEXOS" begin

#         xlsx_in = testpath * "three_nodes.xlsx"
#         xlsx_out = testpath * "three_nodes_PRAS.xlsx"
#         process_workbook(xlsx_in, xlsx_out)

#     end
