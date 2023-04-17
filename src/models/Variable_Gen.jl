const ReEDS_VRE_TYPES = ["wind-ons", "wind-ofs", "dupv", "upv", "distpv", "csp-ns"]

function check_reeds_vre_type(type::Union{STRING, String}) where {STRING <: InlineStrings.InlineString}
    flag = 0
    for vre_type in ReEDS_VRE_TYPES
        if (occursin(vre_type, type))
            flag += 1
        end
    end

    if (flag > 0)
        return true
    else
        return false
    end
end

"""
    This function takes in the attributes of a variable generator (Variable_Gen)
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
    Variable_Gen : Struct
        Returns a struct with all the given attributes.
"""
struct Variable_Gen <: Generator
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
    function Variable_Gen(
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
            error("Check for inconsistencies in $(name) time series data")

        length(capacity) == timesteps ||
            error("The length of the $(name) time series data should be equal
                   to PRAS timesteps")

        check_reeds_vre_type(string(type)) || error("Check the type of $(name) being passed")

        legacy in ["Existing", "New"] || error("Unidentified legacy passed for $(name)")

        0.0 <= FOR <= 1.0 || error("FOR value passed is not allowed for $(name)")

        MTTR > 0 || error("MTTR value passed is not allowed for $(name)")

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

get_capacity(gen::Variable_Gen) = permutedims(round.(Int, gen.capacity))

get_category(gen::Variable_Gen) = "$(gen.legacy)_$(gen.type)"

get_type(gen::Variable_Gen) = gen.type
