"""
Run script for reeds_to_pras routine
"""

include("../src/ReEDS2PRAS.jl")

import PRAS
import ArgParse
import Logging
import LoggingExtras

function parse_commandline()
    """
    Parse command line arguments for use with the reeds_to_pras function
    """
    s = ArgParse.ArgParseSettings()

    @ArgParse.add_arg_table s begin
        "--reedspath"
            help = "Path to ReEDS-2.0 folder"
            required = true
        "--reedscase"
            help = "Path to ReEDS run (usually .../ReEDS-2.0/runs/{casename})"
            required = true
        "--solve_year"
            help = "Year for the case being generated"
            arg_type = Int
            required = true
        "--weather_year"
            help = "The year corresponding to the vg profiles"
            arg_type = Int
            required = true
        "--output_filepath"
            help = "The path for saving the final model. e.g. ./model.pras"
            required = false
        "--timesteps"
            help = "Number of timesteps to use"
            arg_type = Int
            default = 8760
            required = false
        "--log_console"
            help = "Indicate whether to log to console in addition to file"
            arg_type = Int
            default = 1
            required = false
        "--debug"
            help = "Indicate whether to log debug-level messages"
            arg_type = Int
            default = 0
            required = false
    end
    return ArgParse.parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    # Set up logger
    if ~isnothing(parsed_args["output_filepath"])
        logfile = replace(parsed_args["output_filepath"], ".pras"=>".log")

        if parsed_args["debug"] == 1
            logfilehandle = LoggingExtras.MinLevelLogger(
                LoggingExtras.FileLogger(logfile),
                Logging.Debug)
        else
            logfilehandle = LoggingExtras.MinLevelLogger(
                LoggingExtras.FileLogger(logfile),
                Logging.Info)
        end

        if parsed_args["log_console"] == 1
            logger = LoggingExtras.TeeLogger(
                Logging.global_logger(),
                logfilehandle
            )
        else
            logger = logfilehandle
        end

        Logging.global_logger(logger)
    end

    @info "Running reeds_to_pras with the follow inputs:"
    for (arg, val) in parsed_args
        @info "$arg  =>  $val"
    end
    # Run ReEDS2PRAS
    pras_system = ReEDS2PRAS.reeds_to_pras(
        parsed_args["reedscase"], parsed_args["solve_year"],
        parsed_args["reedspath"], parsed_args["timesteps"],
        parsed_args["weather_year"])
    if ~isnothing(parsed_args["output_filepath"])
        PRAS.savemodel(pras_system, parsed_args["output_filepath"])
    end
end

# Run ReEDS2PRAS from command line arguments
main()
