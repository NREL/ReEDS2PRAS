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
function reeds_to_pras(reeds_filepath::String, solve_year::Int64,
                       nems_path::String, timesteps::Int, weather_year::Int)
    """
    Generates a PRAS system from data in ReEDSfilepath

    Parameters
    ----------
    ReEDSfilepath : String
        Location of ReEDS filepath where inputs, results, and outputs are
        stored
    year : Int64
        ReEDS solve year
    NEMS_path : String
        Path to EIA NEMS database
    timesteps : Int
        Number of timesteps
    WEATHERYEAR : Int
        The weather year for variable gen profiles and load

    Returns
    -------
    PRAS.SystemModel
        PRAS SystemModel struct with regions, interfaces, generators,
        region_gen_idxs, storages, region_stor_idxs, generatorstorages,
        region_genstor_idxs, lines, interface_line_idxs, timestamps

    """
    # assume valid weather years as hardcode for now. These should eventually
    # be read in from ReEDS
    if weather_year ∉ [2007, 2008, 2009, 2010, 2011, 2012, 2013]
        msg = "The weather year $weather_year is not a valid VG profile year. "
              "Should be an Int in 2007-2013 currently"
        error(msg)
    end
    ReEDS_data_filepaths = ReEDSdatapaths(reeds_filepath, solve_year)

    @info "creating PRAS system objects..."
    out = load_objects(ReEDS_data_filepaths, weather_year, timesteps,
                       solve_year, nems_path, 2007)
    all_lines, regions, all_gens, storages_array, genstor_array = out
    @info "...objects are created, writing to PRAS system"
    return create_pras_system(regions, all_lines, all_gens, storages_array,
                              genstor_array, timesteps, solve_year)
end


# Run ReEDS2PRAS from command line arguments
main()

end