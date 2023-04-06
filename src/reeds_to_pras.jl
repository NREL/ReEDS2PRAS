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
    reedspath : String
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
function reeds_to_pras(
    reedscase::String,
    solve_year::Int64,
    reedspath::String,
    timesteps::Int,
    weather_year::Int,
)
    # assume valid weather years as hardcode for now. These should eventually
    # be read in from ReEDS
    if weather_year ∉ [2007, 2008, 2009, 2010, 2011, 2012, 2013]
        error(
            "The weather year $weather_year is not a valid VG profile " *
            "year. Should be an Int in 2007-2013 currently",
        )
    end
    ReEDS_data_filepaths = ReEDSdatapaths(reedscase, solve_year)

    @info "creating PRAS system objects..."
    out = load_objects(
        ReEDS_data_filepaths,
        weather_year,
        timesteps,
        solve_year,
        reedspath,
        2007,
    )
    all_lines, regions, all_gens, storages_array, genstor_array = out
    @info "...objects are created, writing to PRAS system"
    return create_pras_system(
        regions,
        all_lines,
        all_gens,
        storages_array,
        genstor_array,
        timesteps,
        weather_year,
    )
end