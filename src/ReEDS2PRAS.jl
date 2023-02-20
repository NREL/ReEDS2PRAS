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
import ArgParse

# Includes
include("load_reeds.jl")
include("process_reeds.jl")
include("create_pras.jl")


function parse_commandline()
    """
    Parse command line arguments for use with the reeds_to_pras function
    """
    s = ArgParse.ArgParseSettings()

    @ArgParse.add_arg_table s begin
        "reeds_filepath"
            help = "Location of ReEDS filepath where inputs, results, and
                    outputs are stored"
            required = true
        "solve_year"
            help = "Year for the case being generated"
            required = true
        "nems_path"
            help = "Path to NEMS public resource database"
            required = true
        "timesteps"
            help = "Number of timesteps to use"
            required = true
        "weather_year"
            help = "The year corresponding to the vg profiles"
            required = true
        "output_filepath"
            help = "The path for saving the final PRAS model. e.g.
                    ./model.pras"
            required = true
    end
    return ArgParse.parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    @info "Running reeds_to_pras with the follow inputs:"
    for (arg, val) in parsed_args
        @info "$arg  =>  $val"
    end
    pras_system = reeds_to_pras(parsed_args["reeds_filepath"],
                                parse(Int64, parsed_args["solve_year"]),
                                parsed_args["nems_path"],
                                parse(Int64, parsed_args["timesteps"]),
                                parse(Int64, parsed_args["weather_year"]))
    PRAS.savemodel(pras_system, parsed_args["output_filepath"])
end


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
        Year for the case being generated
    NEMS_path : String
        Path to NEMS public resource database
    timesteps : Int
        Number of timesteps
    WEATHERYEAR : Int
        The particular year at which vg profiles are based

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