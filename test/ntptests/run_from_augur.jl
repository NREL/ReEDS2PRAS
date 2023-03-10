#%% Inputs
reeds_filepath = "/Users/pbrown/github2/ReEDS-2.0/runs/v20230308_prasM0_USA"
solve_year = 2020
nems_path = "/Users/pbrown/github2/ReEDS-2.0"
timesteps = 8760
weather_year = 2007
output_dir = joinpath(reeds_filepath,"ReEDS_Augur","PRAS")
output_filepath = joinpath(output_dir,"PRAS_m$(solve_year)_w$(weather_year).pras")

#%% Run it
import PRAS
import ArgParse
include("../../src/ReEDS2PRAS.jl")

pras_system = ReEDS2PRAS.reeds_to_pras(
    reeds_filepath, solve_year, nems_path, timesteps, weather_year
)
    
#%% Save it
mkpath(output_dir)
PRAS.savemodel(pras_system, output_filepath)
