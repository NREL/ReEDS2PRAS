# Run Checks for Versioning?
# Could we use this script ot tag our versions because we cannot include ReEDS ver as a dependency?

# Check if you can open a file
function check_file(loc::String)
    io = try
        open(loc)
    catch
        nothing
    end

    if (isnothing(io))
        return nothing, false
    else
        return io, isopen(io)
    end
end

function run_checks(data::ReEDSdatapaths)
    # ReEDS Load Path
    filepath = joinpath(
        data.ReEDSfilepath,
        "ReEDS_Augur",
        "augur_data",
        "pras_load_$(string(data.year)).h5",
    )

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "Load data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the load .h5 file.",
        )
    end

    # Line Capacity Data
    filepath = joinpath(
        data.ReEDSfilepath,
        "ReEDS_Augur",
        "augur_data",
        "tran_cap_$(string(data.year)).csv",
    )

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "Line Capacity data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the tran_cap_$(string(data.year)) .csv file.",
        )
    end

    # Converter Capacity Data
    # Is this Optional?
    filepath = joinpath(
        data.ReEDSfilepath,
        "ReEDS_Augur",
        "augur_data",
        "cap_converter_$(string(data.year)).csv",
    )

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "Converter Capacity data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the cap_converter_$(string(data.year)) .csv file.",
        )
    end

    # Get Technology Types
    # Should we include a check about the structure here? The p1*p10 issue
    filepath = joinpath(data.ReEDSfilepath, "inputs_case", "tech-subset-table.csv")

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "Generator Technology types data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the tech-subset-table .csv file.",
        )
    end

    # Installed Capacity Data
    filepath = joinpath(
        data.ReEDSfilepath,
        "ReEDS_Augur",
        "augur_data",
        "max_cap_$(string(data.year)).csv",
    )

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "Installed Capacity data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the max_cap_$(string(data.year)) .csv file.",
        )
    end

    # Valid Resources
    filepath = joinpath(data.ReEDSfilepath, "inputs_case", "resources.csv")

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "Resource data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the resources .csv file.",
        )
    end

    # Foreced Outage Rate
    filepath = joinpath(
        data.ReEDSfilepath,
        "ReEDS_Augur",
        "augur_data",
        "forced_outage_$(string(data.year)).csv",
    )

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "Forced Outage Rate data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the forced_outage_$(string(data.year)) .csv file.",
        )
    end

    # r-s mapping
    filepath = joinpath(data.ReEDSfilepath, "inputs_case", "rsmap.csv")

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "r-s region mapping data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the rsmap .csv file.",
        )
    end

    # VG CF Data
    filepath = joinpath(
        data.ReEDSfilepath,
        "ReEDS_Augur",
        "augur_data",
        "pras_vre_gen_$(string(data.year)).h5",
    )

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "VG CF data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the pras_vre_gen_$(string(data.year)) .h5 file.",
        )
    end

    # Storage Energy Cap Data
    filepath = joinpath(
        data.ReEDSfilepath,
        "ReEDS_Augur",
        "augur_data",
        "energy_cap_$(string(data.year)).csv",
    )

    io, bool = check_file(filepath)

    if (bool)
        close(io)
    else
        error(
            "Storage Capacity data is not available in ReEDS results. You are either using a ReEDS version not compatible with ReEDS2PRAS (or) the ReEDS case results location passed is erroneous (or) you don't have access to the energy_cap_$(string(data.year)) .csv file.",
        )
    end
end
