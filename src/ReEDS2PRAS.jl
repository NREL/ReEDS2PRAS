"""
Copyright 2022 Alliance for Sustainable Energy and other
NTP Project Developers. See the top-level COPYRIGHT file for details.
Author: Luke Lavin, Surya Chandan Dhulipala
Email: luke.lavin@nrel.gov, suryachandan.dhulipala@nrel.gov
"""
module ReEDS2PRAS
# Exports
export Generators
export Thermal_Gen
export VG_Gen
export Region
export Line

# Imports
import CSV
import DataFrames
import Dates
import HDF5
import PRAS
import Statistics
import TimeZones

# Includes
include("load_reeds.jl")
include("process_reeds.jl")
include("create_pras.jl")

#runs ReEDS2PRAS
function reeds_to_pras(ReEDSfilepath::String, year::Int64, NEMS_path::String, N::Int, WEATHERYEAR::Int)

    if WEATHERYEAR ∉ [2007,2008,2009,2010,2011,2012,2013] # assume valid weather years as hardcode for now. These should eventually be read in from ReEDS
        error("The weather year $WEATHERYEAR is not a valid VG profile year. Should be an Int in 2007-2013 currently")
    end
    ReEDS_data_filepaths = ReEDSdatapaths(ReEDSfilepath,year)
    
    @info "creating PRAS system objects..."
    all_lines,regions,all_gens,storages_array,genstor_array = load_objects(ReEDS_data_filepaths,WEATHERYEAR,N,year,NEMS_path,2007)
    @info "...objects are created, writing to PRAS system"
    return create_pras_system(regions,all_lines,all_gens,storages_array,genstor_array,N,year)
end

end