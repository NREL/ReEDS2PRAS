"Includes functions used for loading ReEDS data"

function clean_names!(input_vec::Vector{<:AbstractString})
    for (idx,a) in enumerate(input_vec)
        if occursin("*",a)
            input_vec[idx] = match(r"\*([a-zA-Z]+-*[a-zA-Z]*)_*", a)[1]
        end
    end
    return input_vec
end

function Load_EIA_NEMS_DB(ReEDS_directory::String)
    EIA_NEMS_loc = joinpath(ReEDS_directory,"inputs","capacitydata","ReEDS_generator_database_final_EIA-NEMS.csv") #there is also a _prm file, not sure which is right?
    EIA_NEMS_data = DataFrames.DataFrame(CSV.File(EIA_NEMS_loc))
    return EIA_NEMS_data
end

struct ReEDSdatapaths
    ReEDSfilepath::String
    year::Int

    # Checks
    function ReEDSdatapaths(x,y)
        (2020 < y <= 2050) || error("year should be between 2020 and 2050 for ReEDS case for now")
        return new(x,y)
    end
end

function get_load_file(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","plot_load_$(string(data.year)).h5")
    isfile(filepath) || error("The year $(data.year) does not have an associated Augur load h5 file. Are you sure ReeDS was run and Augur results saved for $(data.year)?")
    return HDF5.h5read(filepath,"data")
end

function get_vg_cf_data(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","plot_vre_gen_$(string(data.year)).h5")
    isfile(filepath) || error("The year $(data.year) does not have an associated Augur vg h5 file. Are you sure ReeDS was run and Augur results saved for $(data.year)?")
    return HDF5.h5read(filepath,"data")
end

function get_forced_outage_data(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","forced_outage.csv")
    isfile(filepath) || error("No augur-based forced outage data is found.")
    df = DataFrames.DataFrame(CSV.File(filepath,header=true))
    return DataFrames.rename!(df,["ResourceType","FOR"])
end

function get_valid_resources(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"inputs_case","resources.csv")
    isfile(filepath) || error("No resource table is found.")
    return DataFrames.DataFrame(CSV.File(filepath))
end

function get_technology_types(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"inputs_case","tech-subset-table.csv")
    isfile(filepath) || error("no table of technology types!")
    return DataFrames.DataFrame(CSV.File(filepath))
end

function get_line_capacity_data(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","tran_cap_$(string(data.year)).csv") #assumes this file has been formatted by ReEDS to be PRM line capacity data
    isfile(filepath) || error("The year $(data.year) does not have transmission capacity data. Are you sure ReEDS was run and Augur results saved for $(data.year)?")
    return DataFrames.DataFrame(CSV.File(filepath))
end

function get_converter_capacity_data(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","cap_converter_$(string(data.year)).csv")
    isfile(filepath) || error("The year $(data.year) does not have capacity converter data. Are you sure ReEDS was run and Augur results saved for $(data.year)?")
    return DataFrames.DataFrame(CSV.File(filepath))
end

function get_region_mapping(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"inputs_case","rsmap.csv")
    isfile(filepath) || error("no table of r-s region mapping!")
    return DataFrames.DataFrame(CSV.File(filepath))
end

function get_ICAP_data(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","max_cap_$(string(data.year)).csv")
    isfile(filepath) || error("The year $(data.year) does not have generator installed capacity data. Are you sure REEDS was run and Augur results saved for year $(data.year)")
    return DataFrames.DataFrame(CSV.File(filepath))
end

function get_storage_energy_capacity_data(data::ReEDSdatapaths)
    filepath = joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","energy_cap_$(string(data.year)).csv")
    isfile(filepath) || error("The year $(data.year) does not have generator installed storage energy capacity data. Are you sure REEDS was run and Augur results saved for year $(data.year)")
    return DataFrames.DataFrame(CSV.File(filepath))
end