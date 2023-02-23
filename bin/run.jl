"""
Run script for reeds_to_pras routine
"""

include("../src/ReEDS2PRAS.jl")

import PRAS
import ArgParse

function parse_commandline()
    """
    Parse command line arguments for use with the reeds_to_pras function
    """
    s = ArgParse.ArgParseSettings()

    @ArgParse.add_arg_table s begin
        "reeds_filepath"
            help = "Location of ReEDS filepath where inputs, results, and " *
                   "outputs are stored"
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
            help = "The path for saving the final model. e.g. ./model.pras"
            required = false
    end
    return ArgParse.parse_args(s)
end


function main()
    parsed_args = parse_commandline()
    @info "Running reeds_to_pras with the follow inputs:"
    for (arg, val) in parsed_args
        @info "$arg  =>  $val"
    end
    pras_system = ReEDS2PRAS.reeds_to_pras(
        parsed_args["reeds_filepath"], parse(Int64, parsed_args["solve_year"]),
        parsed_args["nems_path"], parse(Int64, parsed_args["timesteps"]),
        parse(Int64, parsed_args["weather_year"]))
    if ~isnothing(parsed_args["output_filepath"])
        PRAS.savemodel(pras_system, parsed_args["output_filepath"])
    end
end

# Run ReEDS2PRAS from command line arguments
main()
