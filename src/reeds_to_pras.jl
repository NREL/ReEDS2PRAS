#runs ReEDS2PRAS
"""
    Generates a PRAS system from data in ReEDSfilepath

    Parameters
    ----------
    ReEDSfilepath : String
        Location of ReEDS filepath where inputs, results, and outputs are
        stored
    year : Int64
        ReEDS solve year
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
function reeds_to_pras(
    reedscase::String,
    solve_year::Int64,
    timesteps::Int,
    weather_year::Int,
)
    # assume valid weather years as hardcode for now. These should eventually
    # be read in from ReEDS
    if weather_year ∉ range(2007,length=7)
        error(
            "The weather year $weather_year is not a valid year for available VG & Load data " *
            "year. Currrently, it should be an Int in [2007-2013].",
        )
    end
    
    ReEDS_data_filepaths = ReEDSdatapaths(reedscase, solve_year)

    @info "Running checks on input data..."
    run_checks(ReEDS_data_filepaths)

    @info "Parsing ReEDS data and creating ReEDS2PRAS objects..."
    out = parse_reeds_data(
        ReEDS_data_filepaths,
        weather_year,
        timesteps,
        solve_year,
        2007,
    )
    lines, regions, gens, storages, genstors = out
    
    @info "ReEDS data successfully parsed, creating a PRAS system"
    return create_pras_system(
        regions,
        lines,
        gens,
        storages,
        genstors,
        timesteps,
        weather_year,
    )
end