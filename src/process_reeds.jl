"""Includes methods and type for processing ReEDS data"""

# Converting FOR and MTTR to λ and μ
"""
    This function calculates the outage rate of a generator based on the
    forecasted generation and its Mean Time To Repair (MTTR).

    Parameters
    ----------
    for_gen: Float64
        fractional amount of units forecasted to be out of service at any given
        time
    mttr: Int64
        Mean Time To Repair of an individual unit

    Returns
    -------
    (λ, μ): Tuple[Float64, Float64]
        λ is the probability of a unit being down, μ is the MTTR in terms of
        number of hours of downtime.
"""
function outage_to_rate(for_gen::Float64, mttr::Int64)
    if (for_gen == 0.0)
        return (λ = 0.0, μ = 1.0)
    else
        if (for_gen > 1.0)
            for_gen = for_gen/100
        end

        if (mttr != 0)
            μ = 1 / mttr
        else
            μ = 0.0
        end
        λ = (μ * for_gen) / (1 - for_gen)

        return (λ = λ, μ = μ)
    end
end


"""
    Constructs the PRAS Region type. Region objects have three main
    attributes - name (String), timesteps (Int64) and load
    (Vector{Float64}). The load attribute represents the region's total
    power demand data over N intervals of measure given in MW, which must
    always be greater than 0.

    Parameters
    ----------
    name : String
        Name to give the region.
    timesteps : Int64
        Number of PRAS timesteps. Must be between 1 and 8784 inclusive.
    load : Vector{Float64}
        Time series data for the region's total power demand must match N
        in length.

    Returns
    -------
    A new instance of the PRAS Region type.
"""
struct Region
    name::String
    timesteps::Int64
    load::Vector{Float64}

    # Inner Constructors & Checks
    function Region(name, timesteps, load=zeros(Float64, timesteps))

        length(load) == timesteps ||
            error("The length of the region load time series data should be
                   equal to PRAS timesteps (timesteps) ")

        all(load .>= 0.0) ||
            error("Check for inconsistencies in load time series data")

        return new(name, timesteps, load)

    end
end

get_name(reg::Region) = reg.name

get_load(reg::Region) = permutedims(round.(Int,reg.load))

abstract type Generator end


"""
    This function is used to define a thermal generator in the model. It
    contains one struct and an inner constructor to check if the inputs are
    valid.

    Parameters
    ----------
    name : String
        The name of the generator.
    timesteps : Int64
        Number of timesteps in the PRAS problem.
    region_name : String
        Name of the region associated with this generator.
    capacity : Float64
        Capacity of the generator.
    fuel : String
        Fuel type of the generator (default "OT").
    legacy : String
        Existing or New generator (default "New").
    FOR : Float64
        Forced Outage Rate (default 0.0).
    MTTR : Int64
        Mean Time To Repair/Replace (default 24).

    Returns
    -------
    An instance of a Thermal_Gen.
"""
struct Thermal_Gen <:Generator
    name::String
    timesteps::Int64
    region_name::String
    capacity::Float64
    fuel::String
    legacy::String
    FOR::Float64
    MTTR::Int64

    # Inner Constructors & Checks
    function Thermal_Gen(name, timesteps, region_name, capacity, fuel="OT",
                         legacy="New", FOR=0.0, MTTR=24)

        capacity >= 0.0 ||
            error("Generator capacity passed is not allowed")

        legacy in ["Existing", "New"] ||
            error("Unidentified legacy passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        return new(name, timesteps, region_name, capacity, fuel, legacy, FOR,
                   MTTR)

    end
end

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
struct VG_Gen <:Generator
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
    function VG_Gen(name, timesteps, region_name, installed_capacity, capacity,
                    type, legacy="New", FOR=0.0, MTTR=24)

        all(0.0 .<= capacity .<= installed_capacity) ||
            error("Check for inconsistencies in VG_Gen time series data")

        length(capacity) == timesteps ||
            error("The length of the VG_Gen time series data should be equal
                   to PRAS timesteps (timesteps) ")

        # type in ["wind-ons","wind-ofs","dupv","upv","csp","distpv"] ||
        #     error("Check the type of VG_Gen being passed")

        legacy in ["Existing","New"] ||
            error("Unidentified legacy passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        return new(name, timesteps, region_name, installed_capacity, capacity, type,
                   legacy, FOR, MTTR)

    end
end

# Getter Functions
get_name(gen::GEN) where {GEN <: Generator} = gen.name

get_legacy(gen::GEN) where {GEN <: Generator} = gen.legacy

get_capacity(gen::Thermal_Gen) = fill(round(Int,gen.capacity),1,gen.timesteps)

get_fuel(gen::Thermal_Gen) = gen.fuel

get_capacity(gen::VG_Gen) = permutedims(round.(Int,gen.capacity))

# Helper Functions
get_outage_rate(gen::GEN) where {GEN <: Generator} = outage_to_rate(gen.FOR,gen.MTTR)

get_λ(gen::GEN) where {GEN <: Generator} = fill(getfield(get_outage_rate(gen),:λ),1,gen.timesteps)

get_μ(gen::GEN) where {GEN <: Generator} = fill(getfield(get_outage_rate(gen),:μ),1,gen.timesteps)

get_category(gen::Thermal_Gen) = "$(gen.legacy)_Thermal_$(gen.fuel)"

get_category(gen::VG_Gen) = "$(gen.legacy)_$(gen.type)"

get_type(gen::VG_Gen) = gen.type

get_type(gen::Thermal_Gen) = gen.fuel

function get_generators_in_region(gens::Vector{<:Generator}, reg_name::String)
    reg_gens = filter(gen -> gen.region_name == reg_name, gens)
    if isempty(reg_gens)
        # @warn "No generators in region: $(reg_name)"
        return Generator[]
    else
        return reg_gens
    end
end

get_generators_in_region(gens::Vector{<:Generator}, reg::Region) = get_generators_in_region(gens, reg.name)

function get_legacy_generators(gens::Vector{<:Generator}, leg::String)
    leg in ["Existing","New"] ||
            error("Unidentified legacy passed")

    leg_gens = filter(gen -> gen.legacy == leg, gens)
    if isempty(leg_gens)
        # @warn "No generators with legacy: $(leg)"
    else
        return leg_gens
    end
end

abstract type Storage end


"""
    This code defines a struct called Battery which is a subtype of
    Storage. The struct has 13 fields: name, timesteps, region_name, type,
    charge_cap, discharge_cap, energy_cap, legacy, charge_eff,
    discharge_eff, carryover_eff, FOR and MTTR. The code also includes an
    inner constructor and checks to ensure that the values passed are
    valid. The constructor checks that timesteps is between 0 and 8784
    (inclusive), charge_cap and discharge_cap are greater than 0.0, energy_cap
    is greater than 0.0, legacy is either "Existing" or "New", all of the
    efficiency values are between 0.0 and 1.0 (inclusive), FOR is between 0.0
    and 1.0 (inclusive) and MTTR is greater than 0. If any of these checks fail
    an error will be thrown.

    Parameters
    ----------
    name : String
        The name of the battery
    timesteps : Int64
        Number of PRAS timesteps
    region_name: String
        Name of the region
    type : String
        Type of battery
    charge_cap : Float64
        Charge capacity
    discharge_cap : Float64
        Discharge capacity
    energy_cap : Float64
        Energy capacity
    legacy : String
        Battery's legacy (existing or new)
    charge_eff : Float64
        Charge efficiency
    discharge_eff : Float64
        Discharge efficiency
    carryover_eff : Float64
        Carryover efficiency
    FOR : Float64
        Factor of restoration
    MTTR : Int64
        Mean time to restore

    Returns
    -------
    Struct with properties related to the batter
"""
struct Battery <:Storage
    name::String
    timesteps::Int64
    region_name::String
    type::String
    charge_cap::Float64
    discharge_cap::Float64
    energy_cap::Float64
    legacy::String
    charge_eff::Float64
    discharge_eff::Float64
    carryover_eff::Float64
    FOR::Float64
    MTTR::Int64

    # Inner Constructors & Checks
    function Battery(name, timesteps, region_name, type, charge_cap,
                     discharge_cap, energy_cap, legacy="New", charge_eff=1.0,
                     discharge_eff=1.0, carryover_eff=1.0, FOR=0.0, MTTR=24)

        @debug "cap_P = $(discharge_cap) MW and cap_E = $(energy_cap) MWh"

        charge_cap > 0.0 ||
            error("Charge capacity passed is not allowed")

        discharge_cap > 0.0 ||
            error("Discharge capacity passed is not allowed")

        energy_cap > 0.0 || error(
            "Energy capacity passed is not allowed: "
            *"$(name) $(charge_cap)MW $(energy_cap)MWh")

        legacy in ["Existing","New"] ||
            error("Unidentified legacy passed")

        all(0.0 .<= [charge_eff, discharge_eff, carryover_eff] .<= 1.0) ||
            error("Invalid charge/discharge/carryover efficiency passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        return new(name, timesteps, region_name, type, charge_cap,
                   discharge_cap, energy_cap, legacy, charge_eff,
                   discharge_eff, carryover_eff, FOR, MTTR)

    end
end


"""
    This code defines a struct called Gen_Storage which is a subtype of
    Storage. The struct has 14 fields: name, timesteps, region_name, type,
    charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap,
    grid_inj_cap, legacy, charge_eff, discharge_eff and carryover_eff. The
    code also contains an inner constructor and checks to ensure that the
    values passed are valid. Specifically, timesteps must be between 0 and 8784
    (inclusive), all capacity values must be greater than or equal to 0.0,
    the legacy value must either be “Existing” or “New”, all of the
    efficiency values must be between 0.0 and 1.0 (inclusive), FOR must be
    between 0.0 and 1.0 (inclusive) and MTTR must be greater than 0.
    Additionally, there is a commented out check that verifies that all of
    the time series data is of the same size. Finally, if any of these
    checks fail, an error will be thrown.

    Parameters
    ----------
    name : string
        Name of Gen_Storage
    timesteps : integer
        PRAS timesteps (timesteps)
    region_name : string
        Region name
    type : string
        Storage type
    charge_cap : float
        Charge capacity
    discharge_cap : float
        Discharge capacity
    energy_cap : float
        Energy capacity
    inflow : float
        Inflow time series data
    grid_withdrawl_cap : float
        Grid withdrawal capacity time series data
    grid_inj_cap : floating point
        Grid injection capacity time series data
    legacy : string
        Must be either "Existing" or "New"
    charge_eff : float
        Charge efficiency
    discharge_eff : float
        Discharge efficiency
    carryover_eff : float
        Carryover efficiency
    FOR : float
        Forced Outage Rate value
    MTTR : integer
        Mean Time To Repair value

    Returns
    -------
    A new instance of Gen_Storage.
"""
struct Gen_Storage <:Storage
    name::String
    timesteps::Int64
    region_name::String
    type::String
    charge_cap::Float64
    discharge_cap::Float64
    energy_cap::Float64
    inflow::Float64
    grid_withdrawl_cap::Float64
    grid_inj_cap::Float64
    legacy::String
    charge_eff::Float64
    discharge_eff::Float64
    carryover_eff::Float64
    FOR::Float64
    MTTR::Int64

    # Inner Constructors & Checks
    function Gen_Storage(name, timesteps, region_name, type, charge_cap,
                         discharge_cap, energy_cap, inflow, grid_withdrawl_cap,
                         grid_inj_cap, legacy ="New", charge_eff=1.0,
                         discharge_eff=1.0, carryover_eff=1.0, FOR=0.0,
                         MTTR=24)

        all(charge_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage charge capacity
                   time series data")

        all(discharge_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage discharge capacity
                   time series data")

        all(energy_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage energy capacity
                   time series data")

        all(inflow .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage inflow time series
                   data")

        all(grid_withdrawl_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage grid withdrawl
                   capacity time series data")

        all(grid_inj_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage grid injection
                   capacity time series data")

        legacy in ["Existing","New"] ||
            error("Unidentified legacy passed")

        all(0.0 .<= [charge_eff, discharge_eff, carryover_eff] .<= 1.0) ||
            error("Invalid charge/discharge/carryover efficiency passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        return new(name, timesteps, region_name, type, charge_cap, discharge_cap,
                   energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap,
                   legacy, charge_eff, discharge_eff, carryover_eff, FOR, MTTR)

    end
end

get_name(stor::STOR) where {STOR <: Storage} = stor.name

get_type(stor::STOR) where {STOR <: Storage} = stor.type

get_legacy(stor::STOR) where {STOR <: Storage} = stor.legacy

get_charge_capacity(stor::Battery) = fill(round(Int, stor.charge_cap), 1, stor.timesteps)

get_charge_capacity(stor::Gen_Storage) = fill(round(Int, stor.charge_cap), 1, stor.timesteps)

get_discharge_capacity(stor::Battery) = fill(round(Int, stor.discharge_cap), 1, stor.timesteps)

get_discharge_capacity(stor::Gen_Storage) = fill(round(Int, stor.discharge_cap), 1, stor.timesteps)

get_energy_capacity(stor::Battery) = fill(round(Int, stor.energy_cap), 1, stor.timesteps)

get_energy_capacity(stor::Gen_Storage) = fill(round(Int ,stor.energy_cap), 1, stor.timesteps)

get_inflow(stor::Gen_Storage) = fill(round(Int, stor.inflow), 1, stor.timesteps)

get_grid_withdrawl_capacity(stor::Gen_Storage) = fill(round(Int, stor.grid_withdrawl_cap), 1, stor.timesteps)

get_grid_injection_capacity(stor::Gen_Storage) = fill(round(Int, stor.grid_inj_cap), 1, stor.timesteps)

get_charge_efficiency(stor::STOR) where {STOR <: Storage} = fill(stor.charge_eff, 1, stor.timesteps)

get_discharge_efficiency(stor::STOR) where {STOR <: Storage} = fill(stor.discharge_eff, 1, stor.timesteps)

get_carryover_efficiency(stor::STOR) where {STOR <: Storage} = fill(stor.carryover_eff, 1, stor.timesteps)

# Helper Functions
get_outage_rate(stor::STOR) where {STOR <: Storage} = outage_to_rate(stor.FOR, stor.MTTR)

get_λ(stor::STOR) where {STOR <: Storage} = fill(getfield(get_outage_rate(stor),:λ), 1, stor.timesteps)

get_μ(stor::STOR) where {STOR <: Storage} = fill(getfield(get_outage_rate(stor),:μ), 1, stor.timesteps)

get_category(stor::STOR) where {STOR <: Storage} = "$(stor.legacy)_$(stor.type)"


"""
    This function searches an array stors of type Vector{<:Storage} for
    storages located in a specific region reg_name. First, it filters the array
    for storages with a region_name field equal to the region name given. If no
    such storages exist, a warning is issued and an empty array of type
    Storage[] is returned. Otherwise, an array containing all the storages from
    this region is returned.

    Parameters
    ----------
    stors : Vector{<:Storage}
        An array of instances of type Storage.
    reg_name : String
        The name of the region to search for in storages.

    Returns
    -------
    reg_stors : Vector{<:Storage}
        An array of Storage instances found in the specified region reg_name.
"""
function get_storages_in_region(stors::Vector{<:Storage}, reg_name::String)
    reg_stors = filter(stor -> stor.region_name == reg_name, stors)
    if isempty(reg_stors)
        # @warn "No storages in region: $(reg_name)"
        return Storage[]
    else
        return reg_stors
    end
end

get_storages_in_region(stors::Vector{<:Storage},
                       reg::Region) = get_storages_in_region(stors, reg.name)


"""
    Get the storage objects which match a given legacy ('Existing' or 'New').

    Parameters
    ----------
    stors: Vector{<:Storage}
        The array of storage objects.
    leg: str
        Legacy of the storage objects. Accepted values are 'Existing' and
        'New'.

    Returns
    -------
    leg_stors: <:Storage
        A subset of ``stors`` that has matching legacy.
        Returns an empty array if there is no match.
"""
function get_legacy_storages(stors, leg)
    leg in ["Existing","New"] ||
            error("Unidentified legacy passed")

    leg_stors = filter(stor -> stor.legacy == leg, stors)
    if isempty(leg_stors)
        # @warn "No storages with legacy: $(leg)"
    else
        return leg_stors
    end
end

emptyvec(::Vector{<:Generator}) = Generator[]

get_components(comps::Vector{<:Generator},
               region_name::String) = get_generators_in_region(comps,
                                                               region_name)

emptyvec(::Vector{<:Storage}) = Storage[]

get_components(comps::Vector{<:Storage},
               region_name::String) = get_storages_in_region(comps,
                                                             region_name)

# Functions for processing ReEDS2PRAS generators and storages
# (preparing PRAS lines)
"""
    Gets components in each region of the system and reorganizes them into
    single sorted component vector and corresponding component index vector.

    Parameters
    ----------
    comps : COMPONENTS
        Vector containing components for each region
    region_names : Vector{String}
        Vector with names of regions in the system

    Returns
    -------
    sorted_comps : COMPONENTS
        Sorted vector of all components from every region
    region_comp_idxs : UnitRange{Int64}, 1
        Index vector pointing to components belonging to each specified region
"""
function get_sorted_components(
        comps::COMPONENTS, region_names::Vector{String}) where
        {COMPONENTS <: Union{Vector{<:Generator}, Vector{<:Storage}}}

    num_regions = length(region_names)
    all_comps = []
    start_idx = Array{Int64}(undef, num_regions);
    region_comp_idxs = Array{UnitRange{Int64}, 1}(undef, num_regions);

    for (idx,region_name) in enumerate(region_names)
        region_comps = get_components(comps, region_name)
        push!(all_comps,region_comps)
        if idx == 1
            start_idx[idx] = 1
        else
            prev_idx = start_idx[idx-1]
            prev_length = length(all_comps[idx-1])
            start_idx[idx] = prev_idx + prev_length
        end
        region_comp_idxs[idx] = range(start_idx[idx],
                                      length=length(all_comps[idx]))
    end

    sorted_comps = emptyvec(comps)
    for idx in eachindex(all_comps)
        if (length(all_comps[idx]) != 0)
            append!(sorted_comps, all_comps[idx])
        end
    end
    return sorted_comps, region_comp_idxs
end

get_sorted_components(
    comps::COMPONENTS, regions::Vector{Region}) where
    {COMPONENTS <: Union{Vector{<:Generator}, Vector{<:Storage}}} =
    get_sorted_components(comps,get_name.(regions))

# Lines
"""
    Constructs a model of LINE.

    Parameters
    ----------
    name : String
        Name chose for a line object
    timesteps : Int64
        Number of timeseries for the same identifier (assumed same for all
        type of objects)
    category : String
        Type of Line object. Category can be one of AC/B2B/LCC/VSC/VSC
        DC-AC converter
    region_from : String
        Origin region of the line object
    region_to : String
        Destination region of the line object
    forward_cap : Float64
        Capacity of the line object in forward direction
    backward_cap : Float64
        Capacity of the line object in backward direction
    legacy : String
        Check whether it is existing or new i.e., RET or TEST
    FOR : Float64
        Forced Outage Rate of the line object
    MTTR : Int64
        Mean Time to Repair of the line object
    VSC : Bool
        Voltage Source Converter of the line object
    converter_capacity : Dict{String, Float64}
        Dictionary that stores the VSC capacities within Region_From and
        Region_To

    Returns
    -------
    Out: Line
        A new instance of Line object as defined above

"""
struct Line
    name::String
    timesteps::Int64
    category::String
    region_from::String
    region_to::String
    forward_cap::Float64
    backward_cap::Float64
    legacy::String
    FOR::Float64
    MTTR::Int64
    VSC::Bool
    converter_capacity:: Dict{String, Float64}

    # Inner Constructors & Checks
    function Line(name, timesteps, category, region_from, region_to,
                  forward_cap, backward_cap, legacy="New", FOR=0.0, MTTR=24,
                  VSC=false, converter_capacity=Dict(region_from => 0.0,
                                                     region_to => 0.0))

        category in ["AC", "B2B", "LCC", "VSC", "VSC DC-AC converter"] ||
            error("Check the category of Line passed")

        ~(region_from == region_to) ||
            error("Region_From and Region_To cannot be the same. PRAS only
                   considers inter-regional lines in Zonal analysis")

        all([forward_cap, backward_cap] .>= 0.0) ||
            error("Check the forward/backward capacity of Line passed")

        legacy in ["Existing","New"] ||
            error("Unidentified legacy passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        all([region_from, region_to] .∈  Ref(keys(converter_capacity))) ||
            error("Check the keys of converter capacity dictionary for VSC DC
                   line")

        return new(name, timesteps, category, region_from, region_to, forward_cap,
                   backward_cap, legacy, FOR, MTTR, VSC, converter_capacity)
    end
end

get_name(ln::Line) = ln.name

get_category(ln::Line) = ln.category

get_forward_capacity(ln::Line) = fill(round(Int,ln.forward_cap),1,ln.timesteps)

get_backward_capacity(ln::Line) = fill(round(Int,ln.backward_cap),1,ln.timesteps)

get_region_from(ln::Line) = ln.region_from

get_region_to(ln::Line) = ln.region_to

#Helper functions
get_outage_rate(ln::Line) = outage_to_rate(ln.FOR,ln.MTTR)

get_λ(ln::Line) = fill(getfield(outage_to_rate(ln.FOR,ln.MTTR),:λ),1,ln.timesteps)

get_μ(ln::Line) = fill(getfield(outage_to_rate(ln.FOR,ln.MTTR),:μ),1,ln.timesteps)


"""
    Return all lines with the specified legacy.

    Parameters
    ----------
    lines : Vector{<:Line}
        List of Lines to filter through.
    leg : String
        Legacy of Line, either 'Existing' or 'New'.

    Returns
    -------
    Vector{<:Line}
        List of all lines with the specified legacy.
"""
function get_legacy_lines(lines::Vector{Line}, leg::String)
    leg in ["Existing","New"] ||
            error("Unidentified legacy passed")

    leg_lines = filter(ln -> ln.legacy == leg, lines)
    if isempty(leg_lines)
        # @warn "No lines with legacy: $(leg)"
    else
        return leg_lines
    end
end

# Functions for  processing ReEDS2PRAS lines (preparing PRAS lines)
"""
    get_legacy_lines function filters the lines Vector{Line} for only those
    lines with a given legacy (either Existing or New) and returns these
    filtered lines.

    Parameters
    ----------
    lines : Vector{Line}
        Vector of lines to be filtered
    leg : String
        the legacy which the lines should have

    Returns
    -------
    leg_lines : Vector{Line}
        Vector of filtered lines with the same legacy as provided
"""
function get_sorted_region_tuples(lines::Vector{Line},
                                  region_names::Vector{String})
    region_idxs = Dict(name => idx for (idx, name) in enumerate(region_names))

    line_from_to_reg_idxs = similar(lines, Tuple{Int,Int})

    for (l, line) in enumerate(lines)

        from_name = get_region_from(line)
        to_name = get_region_to(line)

        from_idx = region_idxs[from_name]
        to_idx = region_idxs[to_name]

        line_from_to_reg_idxs[l] = from_idx < to_idx ? (from_idx, to_idx) :
                                   (to_idx, from_idx)
    end

    return line_from_to_reg_idxs
end

"""
    Returns a list of tuples sorted by region name.

    Parameters
    ----------
    lines : Vector{Line}
        A list of lines containing the regions from and to.

    Returns
    -------
    List[Tuple[str]]
        A list of tuples sorted by region name.
"""
function get_sorted_region_tuples(lines::Vector{Line})
    regions_from = get_region_from.(lines)
    regions_to = get_region_to.(lines)

    region_names = unique(append!(regions_from, regions_to))

    get_sorted_region_tuples(lines, region_names)
end

function get_sorted_region_tuples(lines::Vector{Line}, regions::Vector{Region})
    get_sorted_region_tuples(lines, get_name.(regions))
end

function get_sorted_lines(lines::Vector{Line}, region_names::Vector{String})
    line_from_to_reg_idxs = get_sorted_region_tuples(lines, region_names)
    line_ordering = sortperm(line_from_to_reg_idxs)

    sorted_lines = lines[line_ordering]
    sorted_from_to_reg_idxs = line_from_to_reg_idxs[line_ordering]
    interface_reg_idxs = unique(sorted_from_to_reg_idxs)

    # Ref tells Julia to use interfaces as Vector, only broadcasting over
    # lines_sorted
    interface_line_idxs = searchsorted.(Ref(sorted_from_to_reg_idxs),
                                        interface_reg_idxs)

    return sorted_lines, interface_reg_idxs, interface_line_idxs
end

get_sorted_lines(lines::Vector{Line},
                 regions::Vector{Region}) =
                 get_sorted_lines(lines,get_name.(regions))

"""
    This code takes in a vector of Lines and a vector of Regions as input
    parameters. It filters the Lines to find VSC (voltage source converter)
    lines and non-VSC lines. Then, for each VSC line, it creates two new Line
    objects representing direct current (DC) for the line and two converters
    for the regional connections to the VSC line. Finally, a vector of all
    lines and a vector of all regions are returned as output.

    Parameters
    ----------
    lines : Vector[Line]
        Vector of Line objects that contain information about all the lines in
        the system
    regions : Vector[Region]
        Vector of Region objects that contain information about all regions in
        the system

    Returns
    -------
    non_vsc_dc_lines : Vector[Line]
        Vector of Line objects with non_vsc_dc_lines
    regions : Vector[Region]
        Vector of Region objects with added DC lines
"""
function process_vsc_lines(lines::Vector{Line}, regions::Vector{Region})
    timesteps = first(regions).timesteps
    non_vsc_dc_lines = filter(line -> ~line.VSC, lines)
    vsc_dc_lines = filter(line -> line.VSC, lines)

    for vsc_line in vsc_dc_lines
        dc_region_from = "DC_$(vsc_line.region_from)"
        dc_region_to = "DC_$(vsc_line.region_to)"

        for reg_name in [dc_region_from, dc_region_to]
            if ~(reg_name in get_name.(regions))
                push!(regions, Region(reg_name, timesteps,
                                      zeros(Float64, timesteps)))
            end
        end

        push!(non_vsc_dc_lines,
              Line("$(vsc_line.name)_DC", vsc_line.timesteps,
                   vsc_line.category, dc_region_from, dc_region_to,
                   vsc_line.forward_cap, vsc_line.backward_cap,
                   vsc_line.legacy, vsc_line.FOR, vsc_line.MTTR))
        push!(non_vsc_dc_lines,
              Line("$(vsc_line.name)_Converter_From", vsc_line.timesteps,
                   vsc_line.category, dc_region_from, vsc_line.region_from,
                   vsc_line.converter_capacity[vsc_line.region_from],
                   vsc_line.converter_capacity[vsc_line.region_from],
                   vsc_line.legacy,vsc_line.FOR, vsc_line.MTTR))
        push!(non_vsc_dc_lines,
              Line("$(vsc_line.name)_Converter_To", vsc_line.timesteps,
                   vsc_line.category, dc_region_to, vsc_line.region_to,
                   vsc_line.converter_capacity[vsc_line.region_to],
                   vsc_line.converter_capacity[vsc_line.region_to],
                   vsc_line.legacy,vsc_line.FOR, vsc_line.MTTR))

    end
    return non_vsc_dc_lines, regions
end