
#need to figure out what the appropriate way is to load in the right kind of .h5 file to get this to work...

function process_solution(ReEDSfilepath::String,inputpath_h5::String)
    #,timestep::Period=Hour(1),timezone::TimeZone=tz"UTC"
    h5open(inputpath_h5, "r") do prasfile::HDF5.File
        # process_metadata!(prasfile,timestep,timezone)
        #, timestep, timezone
        process_regions!(ReEDSfilepath,prasfile)

    end

end

function process_metadata!(
    prasfile::HDF5.File,
    timestep::Period,
    timezone::TimeZone)

    # version_message = "Only H5PLEXOS v0.6 files are supported"
    # if haskey(attributes(plexosfile), "h5plexos")
    #     version = read(attributes(plexosfile)["h5plexos"])
    #     version_match = match(r"^v0.6.\d+$", version)
    #     isnothing(version_match) && error(version_message * ", got " * version)
    # else
    #     error(version_message)
    # end

    attrs = attributes(prasfile)
end

function process_regions!(ReEDSfilepath::String,prasfile::HDF5.File)
    #stringlength::Int, compressionlevel::Int

    # Load required data from plexosfile
    regiondata = HDF5.h5read(joinpath(ReEDSfilepath,"inputs_case","load.h5"), "data") 
    # load_data = 
    # regions = 
    regions = regiondata["axis0"]
    load = regiondata["block0_values"]

    n_regions = length(regiondata)

    # Save data to prasfile
    regions = create_group(prasfile, "regions")
    # string_table!(regions, "_core", regiondata[!, [:name]], stringlength)
    # regions["load", compress=compressionlevel] = round.(UInt32, load)

    # return n_regions
end