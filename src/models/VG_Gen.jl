"""
    This function takes in the attributes of a variable generator (VG_Gen)
    and returns an object containing all its information. The input
    parameters are:

    Parameters
    ----------
    name : String
        Name of the Variable Generator.
    timesteps : Int64
        Number of timesteps for the PRAS model.
    region_name : String
        Name of the Region where Variable Generator's load area is present.
    installed_capacity : Float64
        Installed Capacity of the Variable Generator.
    capacity : Vector{Float64}
        Capacity factor time series data ('forecasted capacity' /
        'nominal capacity') for the PRAS Model.
    type : String
        Type of Variable Generator being passed.
    legacy : String
        State of the Variable Generator, i.e., existing or new.
    FOR : Float64
        Forced Outage Rate parameter of the Variable Generator.
    MTTR : Int64
        Mean Time To Repair parameter of the Variable Generator.

    Returns
    -------
    VG_Gen : Struct
        Returns a struct with all the given attributes.
"""
struct VG_Gen <: Generator
    name::String
    timesteps::Int64
    region_name::String
    installed_capacity::Float64
    capacity::Vector{Float64}
    type::String
    legacy::String
    FOR::Float64
    MTTR::Int64

    # Inner Constructors & Checks
    function VG_Gen(
        name,
        timesteps,
        region_name,
        installed_capacity,
        capacity,
        type,
        legacy = "New",
        FOR = 0.0,
        MTTR = 24,
    )
        all(0.0 .<= capacity .<= installed_capacity) ||
            error("Check for inconsistencies in VG_Gen time series data")

        length(capacity) == timesteps ||
            error("The length of the VG_Gen time series data should be equal
                   to PRAS timesteps (timesteps) ")

        # type in ["wind-ons","wind-ofs","dupv","upv","csp","distpv"] ||
        #     error("Check the type of VG_Gen being passed")

        legacy in ["Existing", "New"] || error("Unidentified legacy passed")

        0.0 <= FOR <= 1.0 || error("FOR value passed is not allowed")

        MTTR > 0 || error("MTTR value passed is not allowed")

        return new(
            name,
            timesteps,
            region_name,
            installed_capacity,
            capacity,
            type,
            legacy,
            FOR,
            MTTR,
        )
    end
end

# Getter Functions

get_capacity(gen::VG_Gen) = permutedims(round.(Int, gen.capacity))

get_category(gen::VG_Gen) = "$(gen.legacy)_$(gen.type)"

get_type(gen::VG_Gen) = gen.type
