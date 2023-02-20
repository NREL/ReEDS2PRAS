"Includes functions used for loading ReEDS data"

function clean_names!(input_vec::Vector{<:AbstractString})
    """
    Loop through the vector and replace any '*' present in the elements with
    the text between the asterisk, excluding the '_'.

    Parameters
    ----------
    input_vec : Vector{<:AbstractString}
        Vector of strings that need to be parsed

    Returns
    -------
    input_vec : Vector{<:AbstractString}
        Vector of strings which has been cleaned of '*'
    """
    for (idx,a) in enumerate(input_vec)
        if occursin("*",a)
            input_vec[idx] = match(r"\*([a-zA-Z]+-*[a-zA-Z]*)_*", a)[1]
        end
    end
    return input_vec
end


function Load_EIA_NEMS_DB(ReEDS_directory::String)
    """
    Loads the EIA-NEMS Generator Database from a given ReEDS directory.

    Parameters
    ----------
    ReEDS_directory : String
        Path of the ReEDS directory from which to load the EIA-NEMS Generator
        Database data.

    Returns
    -------
    EIA_NEMS_data : DataFrame
        A DataFrame containing the EIA-NEMS Generator Database data.
    """
    #there is also a _prm file, not sure which is right?
    EIA_NEMS_loc = joinpath(ReEDS_directory, "inputs", "capacitydata",
                            "ReEDS_generator_database_final_EIA-NEMS.csv")
    EIA_NEMS_data = DataFrames.DataFrame(CSV.File(EIA_NEMS_loc))
    return EIA_NEMS_data
end


struct ReEDSdatapaths
    ReEDSfilepath::String  # The filepath where Augur data is saved
    year::Int              # year 2020-2050

    function ReEDSdatapaths(x,y)
        """
        Creates a datapath for a given year within the valid date period: 2020
        < y <= 2050.  Used as a parameter for other functions in order to
        access correctly dated input files.

        Parameters
        ----------
        x : String
            Filepath where Augur data is saved
        y : Int
            Year of the case being run

        Returns
        -------
        A new object with filepath and valid year parameters

        """
        msg = ("year should be between 2020 and 2050 for ReEDS case for now")
        (2020 < y <= 2050) || error(msg)
        return new(x,y)
    end
end

function get_load_file(data::ReEDSdatapaths)
    """
    Loads an .h5 file from the given file path, containing Electrical Demand
    data from the indicated year.

    Parameters
    ----------
    data : ReEDSdatapaths
        An instance of ReEDSdatapaths with the necessary arguments set.

    Returns
    -------
    HDF5.h5read(filepath, "data")
        A readout of the Augur load h5 file associated with the given ReEDS
        filepath and year.

    """
    filepath = joinpath(data.ReEDSfilepath, "ReEDS_Augur", "augur_data",
                        "plot_load_$(string(data.year)).h5")
    msg = "The year $(data.year) does not have an associated Augur load h5 "
          "file. Are you sure ReeDS was run and Augur results saved for "
          "$(data.year)?"
    isfile(filepath) || error(msg)
    return HDF5.h5read(filepath,"data")
end


function get_vg_cf_data(data::ReEDSdatapaths)
    """
    This function reads a hdf5 file from the ReEDS Augur directory, based on
    the year provided in the ReEDSdatapaths struct.

    Parameters
    ----------
    data : ReEDSdatapaths
        Struct containing the `ReEDSfilepath` and a year, indicating which .h5
        file should be read.

    Returns
    -------
    The requested ``hdf5`` as a data frame.

    Raises
    ------
    error
        If the year does not have an associated Augur vg h5 file.
    """
    filepath = joinpath(data.ReEDSfilepath, "ReEDS_Augur", "augur_data",
                        "plot_vre_gen_$(string(data.year)).h5")
    msg = "The year $(data.year) does not have an associated Augur vg h5 file.
           Are you sure ReeDS was run and Augur results saved for
           $(data.year)?"
    isfile(filepath) || error(msg)
    return HDF5.h5read(filepath,"data")
end


function get_forced_outage_data(data::ReEDSdatapaths)
    """
    Get the forced outage data from the augur files.

    Parameters
    ----------
    data : ReEDSdatapaths
        Struct containing relevant datapaths and year from which to extract
        the data.

    Returns
    -------
    DataFrames.DataFrame
        Dataframe containing the forced outage data.
    """

    filepath = joinpath(data.ReEDSfilepath, "ReEDS_Augur", "augur_data",
                        "forced_outage.csv")
    isfile(filepath) || error("No augur-based forced outage data is found.")
    df = DataFrames.DataFrame(CSV.File(filepath,header=true))
    return DataFrames.rename!(df,["ResourceType","FOR"])
end


function get_valid_resources(data::ReEDSdatapaths)
    """
    Get the valid resources from a ReEDS case.

    Parameters
    ----------
    data : ReEDSdatapaths
        Struct containing ReEDS filepaths and year

    Returns
    -------
    DataFrames.DataFrame
        A DataFrame containing the valid resources of the ReEDS case
    """
    filepath = joinpath(data.ReEDSfilepath, "inputs_case", "resources.csv")
    isfile(filepath) || error("No resource table is found.")
    return DataFrames.DataFrame(CSV.File(filepath))
end


function get_technology_types(data::ReEDSdatapaths)
    """
    This function gets the technology types from a csv file.

    Arguments
    ---------
    data : ReEDSdatapaths
        A struct containing paths and dates related to ReEDS analyses.

    Returns
    -------
    DataFrames.DataFrame
        The technology type table in the form of a DataFrames object.
    """
    filepath = joinpath(data.ReEDSfilepath, "inputs_case",
                        "tech-subset-table.csv")
    isfile(filepath) || error("no table of technology types!")
    return DataFrames.DataFrame(CSV.File(filepath))
end


function get_line_capacity_data(data::ReEDSdatapaths)
    """
    Gets line capacity data for the given ReEDS database.

    Parameters
    ----------
    data : ReEDSdatapaths
        Contains the filepath of data and year of analysis

    Returns
    -------
    DataFrame
        A dataframe with transmission capacity data; assumes this file has been
        formatted by ReEDS
    """
    #assumes this file has been formatted by ReEDS to be PRM line capacity data
    filepath = joinpath(data.ReEDSfilepath, "ReEDS_Augur", "augur_data",
                        "tran_cap_$(string(data.year)).csv")
    msg = "The year $(data.year) does not have transmission capacity data.
           Are you sure ReEDS was run and Augur results saved for
           $(data.year)?"
    isfile(filepath) || error(msg)
    return DataFrames.DataFrame(CSV.File(filepath))
end


function get_converter_capacity_data(data::ReEDSdatapaths)
    """
    Get the converter capacity data associated with the given ReEDSdatapaths
    object.

    Parameters
    ----------
    data : ReEDSdatapaths)
        A ReEDSdatapaths object containing the relevant file paths and year.

    Returns
    -------
    DataFrames.DataFrame
        The DataFrame of the converter capacity data for the given year.
    """
    filepath = joinpath(data.ReEDSfilepath, "ReEDS_Augur", "augur_data",
                        "cap_converter_$(string(data.year)).csv")
    msg = "The year $(data.year) does not have capacity converter data. Are
           you sure ReEDS was run and Augur results saved for $(data.year)?"
    isfile(filepath) || error(msg)
    return DataFrames.DataFrame(CSV.File(filepath))
end


function get_region_mapping(data::ReEDSdatapaths)
    """
    Returns a DataFrame containing the r-s region mapping from the
    ReEDSdatapaths object.

    Parameters
    ----------
    data : ReEDSdatapaths
        An object containing the filepaths to the ReEDS input files.

    Returns
    -------
    DataFrame
        A DataFrame containing the r-s region mapping.

    Raises
    ------
    Error
        If no table of r-s region mapping is found.

    """
    filepath = joinpath(data.ReEDSfilepath,"inputs_case","rsmap.csv")
    isfile(filepath) || error("no table of r-s region mapping!")
    return DataFrames.DataFrame(CSV.File(filepath))
end


function get_ICAP_data(data::ReEDSdatapaths)
    """
    Returns a DataFrame containing the installed capacity of generators for a
    given year.

    Parameters
    ----------
    data : ReEDSdatapaths
        A ReEDSdatapaths object containing the year and filepath.

    Returns
    -------
    DataFrame
        A DataFrame containing the installed capacity data.

    Raises
    ------
    Error
        If the year does not have generator installed capacity data.
    """
    filepath = joinpath(data.ReEDSfilepath, "ReEDS_Augur", "augur_data",
                        "max_cap_$(string(data.year)).csv")
    msg = "The year $(data.year) does not have generator installed capacity
           data. Are you sure REEDS was run and Augur results saved for year
           $(data.year)"
    isfile(filepath) || error(msg)
    return DataFrames.DataFrame(CSV.File(filepath))
end


function get_storage_energy_capacity_data(data::ReEDSdatapaths)
    """
    Returns a DataFrame containing the installed storage energy capacity data
    for the year specified in the ReEDSdatapaths object.

    Parameters
    ----------
    data : ReEDSdatapaths
        An object containing paths to ReEDS data.

    Returns
    -------
    DataFrame
        DataFrame containing storage energy capacity data

    Raises
    ------
    Error
        If the filepath for the specified year does not exist.
    """
    filepath = joinpath(data.ReEDSfilepath, "ReEDS_Augur", "augur_data",
                        "energy_cap_$(string(data.year)).csv")
    msg = "The year $(data.year) does not have generator installed storage
           energy capacity data. Are you sure REEDS was run and Augur
           results saved for year $(data.year)"
    isfile(filepath) || error(msg)
    return DataFrames.DataFrame(CSV.File(filepath))
end