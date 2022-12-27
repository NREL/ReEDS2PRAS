#######################################################
"""
Copyright 2022 Alliance for Sustainable Energy and other
NTP Project Developers. See the top-level COPYRIGHT file for details.
Author: Surya Chandan Dhulipala, Luke Lavin
Email: suryachandan.dhulipala@nrel.gov, luke.lavin@nrel.gov
"""
# November 2022
#Building a PRAS model based on input ReEDS case 
# Do we need EIA-860 plant data for PV & Wind generators and other ReEDS data
#######################################################
module ReEDS2PRAS
#################################################################################
# Exports
#################################################################################
export generators
export thermal_gen
export vg_gen
export region
export line
# export  make_pras_system_from_mapping_info
#################################################################################
# Imports
#################################################################################
# import HTTP
# import JSON
# import ArchGDAL
import DataFrames
# import XLSX
import CSV
# import Random
import PRAS
import Dates
import HDF5
import Statistics
import TimeZones
#################################################################################
# Includes
#################################################################################
# include("main/naerm_data_service.jl") 
include("main/utils.jl")
include("main/types.jl")
include("main/pras_system_from_mapping.jl")

end