#####################################################
# Surya
# NREL
# December 2022
# ReEDS2PRAS - NTPS
# Structs to make objects from ReEDS data
#######################################################
# Converting FOR and MTTR to λ and μ
#######################################################
function outage_to_rate(for_gen::Float64,mttr::Int64)
    if (for_gen == 0.0)
        return (λ = 0.0, μ = 1.0)
    else
        if (for_gen >1.0)
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
# Regions
struct Region
    name::String
    N::Int64
    load::Vector{Float64}

    # Inner Constructors
    # For illustration purposes
    Region(name, N) = Region(name,N,zeros(Float64,N))

    Region(nothing) = Region("region_1",10,fill(100.0,10))

    # Checks

    function Region(name, N, load)

        0 < N <= 8784 ||
            error("Check the PRAS timesteps (N) passed.")

        length(load) == N ||
            error("The length of the region load time series data should be equal to PRAS timesteps (N)")

        all(load .>= 0.0) ||
            error("Check for inconsistencies in load time series data")

        return new(name, N, load)

    end
end

# Getter Functions (to make broadcasting easier??)
# According to Gord's comment - Do we need these for Region?

get_name(reg::Region) = reg.name

get_load(reg::Region) =  permutedims(round.(Int,reg.load))

# Generators
abstract type Generator end
const Generators = Vector{<:Generator}

struct Thermal_Gen <:Generator
    name::String 
    N::Int64
    region_name::String
    capacity::Float64
    fuel::String
    legacy::String
    FOR::Float64
    MTTR::Int64

    # Inner Constructors
    Thermal_Gen(name, N, region_name, capacity, fuel, legacy, FOR) = Thermal_Gen(name, N, region_name, capacity, fuel, legacy, FOR, 24)

    Thermal_Gen(name, N, region_name, capacity, fuel, legacy) = Thermal_Gen(name, N, region_name, capacity, fuel, legacy, 0.0, 24)

    Thermal_Gen(name, N, region_name, capacity, fuel) = Thermal_Gen(name, N, region_name, capacity, fuel, "New", 0.0, 24)

    Thermal_Gen(name, N, region_name, capacity) = Thermal_Gen(name, N, region_name, capacity, "OT", "New", 0.0, 24)

    # For illustration purposes
    Thermal_Gen(nothing) = Thermal_Gen("thermal_gen_1",10,"reg_1",10.0,"NG","New",0.1,24)

    # Checks

    function Thermal_Gen(name, N, region_name, capacity, fuel, legacy, FOR, MTTR)

        0 < N <= 8784 ||
            error("Check the PRAS timesteps (N) passed.")

        capacity >= 0.0 ||
            error("Generator capacity passed is not allowed")

        legacy in ["Existing","New"] ||
            error("Unidentified legacy passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        return new(name, N, region_name, capacity, fuel, legacy, FOR, MTTR)

    end
end

struct VG_Gen <:Generator
    name::String 
    N::Int64
    region_name::String
    installed_capacity::Float64
    capacity::Vector{Float64}
    type::String
    legacy::String
    FOR::Float64
    MTTR::Int64

    # Inner Constructors
    VG_Gen(name, N, region_name, installed_capacity, capacity, type, legacy) = VG_Gen(name, N, region_name, installed_capacity, capacity, type, legacy, 0.0, 24)
    VG_Gen(name, N, region_name, installed_capacity, capacity, type) =  VG_Gen(name, N, region_name, installed_capacity, capacity, type, "New", 0.0, 24)
   
    # Checks
    function VG_Gen(name, N, region_name, installed_capacity, capacity, type, legacy, FOR, MTTR)

        0 < N <= 8784 ||
            error("Check the PRAS timesteps (N) passed.")

        all(0.0 .<= capacity .<= installed_capacity) ||
            error("Check for inconsistencies in VG_Gen time series data")

        length(capacity) == N ||
            error("The length of the VG_Gen time series data should be equal to PRAS timesteps (N)")

        # type in ["wind-ons","wind-ofs","dupv","upv","csp","distpv"] ||
        #     error("Check the type of VG_Gen being passed")

        legacy in ["Existing","New"] ||
            error("Unidentified legacy passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        return new(name, N, region_name, installed_capacity, capacity, type, legacy, FOR, MTTR)

    end
end

# Getter Functions
get_name(gen::GEN) where {GEN <: Generator} = gen.name

get_legacy(gen::GEN) where {GEN <: Generator} = gen.legacy

get_capacity(gen::Thermal_Gen) = fill(round(Int,gen.capacity),1,gen.N)

get_fuel(gen::Thermal_Gen) = gen.fuel

get_capacity(gen::VG_Gen) = permutedims(round.(Int,gen.capacity))

# Helper Functions
get_outage_rate(gen::GEN) where {GEN <: Generator} = outage_to_rate(gen.FOR,gen.MTTR)

get_λ(gen::GEN) where {GEN <: Generator} = fill(getfield(get_outage_rate(gen),:λ),1,gen.N)

get_μ(gen::GEN) where {GEN <: Generator} = fill(getfield(get_outage_rate(gen),:μ),1,gen.N)

get_category(gen::Thermal_Gen) = "$(gen.legacy)_Thermal_$(gen.fuel)"

get_type(gen::VG_Gen) = gen.type

get_type(gen::Thermal_Gen) = gen.fuel

get_category(gen::VG_Gen) = "$(gen.legacy)_$(gen.type)"

function get_generators_in_region(gens::Generators, reg_name::String)
    reg_gen_idxs = findall(getfield.(gens,:region_name) .== reg_name)
    if isnothing(reg_gen_idxs)
        @warn "No generators in the region"
    else
        return gens[reg_gen_idxs]
    end
end

get_generators_in_region(gens::Generators, reg::Region) = get_generators_in_region(gens, reg.name)

function get_legacy_generators(gens::Generators, leg::String)
    leg in ["Existing","New"] || error("Unidentified legacy passed")

    leg_gen_idxs = findall(getfield.(gens,:legacy) .== leg)
    if isnothing(leg_gen_idxs)
        @warn "No generators with this legacy"
    else
        return gens[leg_gen_idxs]
    end
end

# Storages
abstract type Storage end
const Storages = Vector{<:Storage}

struct Battery <:Storage
    name::String 
    N::Int64
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

    # Inner Constructors
    Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, legacy, charge_eff, discharge_eff, carryover_eff) = 
                                       Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, legacy, charge_eff, discharge_eff, carryover_eff, 0.0, 24)

    Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, legacy, charge_eff, discharge_eff) = 
                                       Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, legacy, charge_eff, discharge_eff, 1.0, 0.0, 24)

    Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, legacy) = 
                                       Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, legacy, 1.0, 1.0, 1.0, 0.0, 24)

    Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap) = 
                                       Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, "New", 1.0, 1.0, 1.0, 0.0, 24)
    
    # Checks
    function Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, legacy, charge_eff, discharge_eff, carryover_eff, FOR, MTTR)

        0 < N <= 8784 ||
            error("Check the PRAS timesteps (N) passed.")

        charge_cap > 0.0 ||
            error("Charge capacity passed is not allowed")

        discharge_cap > 0.0 ||
            error("Discharge capacity passed is not allowed")

        energy_cap > 0.0 ||
            error("Energy capacity passed is not allowed")

        legacy in ["Existing","New"] ||
            error("Unidentified legacy passed")

        all(0.0 .<= [charge_eff, discharge_eff, carryover_eff] .<= 1.0) ||
            error("Invalid charge/discharge/carryover efficiency passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        return new(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, legacy, charge_eff, discharge_eff, carryover_eff, FOR, MTTR)
    end
end

struct Gen_Storage <:Storage
    name::String 
    N::Int64
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

    # Inner Constructors
    Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, legacy, charge_eff, discharge_eff, 
                carryover_eff) = Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, legacy, charge_eff, discharge_eff, 
                                             carryover_eff, 0.0, 24)

    Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, legacy, charge_eff, discharge_eff) = 
                Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, legacy, charge_eff, 
                            discharge_eff, 1.0, 0.0, 24)

    Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, legacy) = 
                Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, legacy, 1.0, 1.0, 1.0, 0.0, 24)

    Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap) = 
                Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, "New", 1.0, 1.0, 1.0, 0.0, 24)
    
    # Checks
    function Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, legacy, charge_eff, discharge_eff, 
                         carryover_eff, FOR, MTTR)

        0 < N <= 8784 ||
            error("Check the PRAS timesteps (N) passed.")

        all(charge_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage charge capacity time series data")

        all(discharge_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage discharge capacity time series data")

        all(energy_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage energy capacity time series data")

        all(inflow .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage inflow time series data")

        all(grid_withdrawl_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage grid withdrawl capacity time series data")

        all(grid_inj_cap .>= 0.0) ||
            error("Check for inconsistencies in Gen_Storage grid injection capacity time series data")

        # all(length.([charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap]) .== N) ||
        #     error("Invalid charge/discharge/carryover efficiency passed")
            
        legacy in ["Existing","New"] ||
            error("Unidentified legacy passed")

        all(0.0 .<= [charge_eff, discharge_eff, carryover_eff] .<= 1.0) ||
            error("Invalid charge/discharge/carryover efficiency passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        return new(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, legacy, charge_eff, discharge_eff, 
                   carryover_eff, FOR, MTTR)

    end
end

get_name(stor::STOR) where {STOR <: Storage} = stor.name

get_type(stor::STOR) where {STOR <: Storage} = stor.type

get_legacy(stor::STOR) where {STOR <: Storage} = stor.legacy

get_outage_rate(stor::STOR) where {STOR <: Storage} = outage_to_rate(stor.FOR,stor.MTTR)

get_λ(stor::STOR) where {STOR <: Storage} = fill(getfield(outage_to_rate(stor.FOR,stor.MTTR),:λ),1,stor.N)

get_μ(stor::STOR) where {STOR <: Storage} = fill(getfield(outage_to_rate(stor.FOR,stor.MTTR),:μ),1,stor.N)

get_category(stor::STOR) where {STOR <: Storage} = "$(stor.legacy)_$(stor.type)"

get_charge_capacity(stor::Battery) = fill(round(Int,stor.charge_cap),1,stor.N)

get_charge_capacity(stor::Gen_Storage) = fill(round(Int,stor.charge_cap),1,stor.N)#permutedims(round.(Int,stor.charge_cap))

get_discharge_capacity(stor::Battery) = fill(round(Int,stor.discharge_cap),1,stor.N)

get_discharge_capacity(stor::Gen_Storage) = fill(round(Int,stor.discharge_cap),1,stor.N)#permutedims(round.(Int,stor.discharge_cap))

get_energy_capacity(stor::Battery) = fill(round(Int,stor.energy_cap),1,stor.N)

get_energy_capacity(stor::Gen_Storage) = fill(round(Int,stor.energy_cap),1,stor.N)#permutedims(round.(Int,stor.energy_cap))

get_inflow(stor::Gen_Storage) = fill(round(Int,stor.inflow),1,stor.N)#permutedims(round.(Int,stor.inflow))

get_grid_withdrawl_capacity(stor::Gen_Storage) = fill(round(Int,stor.grid_withdrawl_cap),1,stor.N)#permutedims(round.(Int,stor.grid_withdrawl_cap))

get_grid_injection_capacity(stor::Gen_Storage) = fill(round(Int,stor.grid_inj_cap),1,stor.N)#permutedims(round.(Int,stor.grid_inj_cap))

get_charge_efficiency(stor::STOR) where {STOR <: Storage} = fill(stor.charge_eff,1,stor.N)

get_discharge_efficiency(stor::STOR) where {STOR <: Storage} = fill(stor.discharge_eff,1,stor.N)

get_carryover_efficiency(stor::STOR) where {STOR <: Storage} = fill(stor.carryover_eff,1,stor.N)

function get_storages_in_region(stors::Storages, reg_name::String)
    reg_stor_idxs = findall(getfield.(stors,:region_name) .== reg_name)
    if isnothing(reg_stor_idxs)
        @warn "No storages in the region"
    else
        return stors[reg_stor_idxs]
    end
end

get_storages_in_region(stors::Storages, reg::Region) = get_storages_in_region(stors, reg.name)

function get_legacy_storages(stors::Storages, leg::String)
    leg in ["Existing","New"] || error("Unidentified legacy passed")

    leg_stor_idxs = findall(getfield.(stors,:legacy) .== leg)
    if isnothing(eg_stor_idxs)
        @warn "No storages with this legacy"
    else
        return stors[leg_stor_idxs]
    end
end

component_dict = Dict(
    Vector{Generator} => (func = get_generators_in_region, vec = Generator[]),
    Vector{Storage} => (func = get_storages_in_region, vec = Storage[]),
    Vector{Gen_Storage} => (func = get_storages_in_region, vec = Storage[])
)
# Functions for processing ReEDS2PRAS generators and storages (preparing PRAS lines)
function get_sorted_components(comps::COMPONENTS, region_names::Vector{String}) where {COMPONENTS <: Union{Generators,Storages}}
    num_regions = length(region_names)
    all_comps = []
    start_id = Array{Int64}(undef,num_regions); 
    region_comp_idxs = Array{UnitRange{Int64},1}(undef,num_regions); 

    for (idx,region_name) in enumerate(region_names)
        region_comps = getfield(component_dict[typeof(comps)],:func)(comps,region_name)
        push!(all_comps,region_comps)
        idx==1 ? start_id[idx] = 1 : start_id[idx] =start_id[idx-1]+length(all_comps[idx-1])
        region_comp_idxs[idx] = range(start_id[idx], length=length(all_comps[idx]))
    end
    
    sorted_comps = getfield(component_dict[typeof(comps)],:vec)
    for idx in eachindex(all_comps)
        if (length(all_comps[idx]) != 0)
            append!(sorted_comps,all_comps[idx])
        end
    end
    return sorted_comps, region_comp_idxs
end

get_sorted_components(comps::COMPONENTS, regions::Vector{Region}) where {COMPONENTS <: Union{Generators,Storages}} = get_sorted_components(comps,get_name.(regions))

# Lines
struct Line
    name::String 
    N::Int64
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
    
    # Inner Constructors
    Line(name, N, category, region_from, region_to, forward_cap, backward_cap, legacy, FOR, MTTR) = 
                            Line(name, N, category, region_from, region_to, forward_cap, backward_cap, legacy, FOR, MTTR, false, Dict(region_from => 0.0, region_to => 0.0))

    Line(name, N, category, region_from, region_to, forward_cap, backward_cap, legacy) = 
                            Line(name, N, category, region_from, region_to, forward_cap, backward_cap, legacy, false, Dict(region_from => 0.0, region_to => 0.0))

    Line(name, N, category, region_from, region_to, forward_cap, backward_cap) = 
                            Line(name, N, category, region_from, region_to, forward_cap, backward_cap, "New", 0.0, 24, false,  Dict(region_from => 0.0, region_to => 0.0))

    Line(name, N, category, region_from, region_to, forward_cap) =
                            Line(name, N, category, region_from, region_to, forward_cap, forward_cap, "New", 0.0, 24, false, Dict(region_from => 0.0, region_to => 0.0))

    # Checks
    function Line(name, N, category, region_from, region_to, forward_cap, backward_cap, legacy, FOR, MTTR, VSC, converter_capacity)

        0 < N <= 8784 ||
            error("Check the PRAS timesteps (N) passed.")

        category in ["AC","Two-Terminal","VSC","VSC DC-AC converter"] ||
            error("Check the category of Line passed")

        ~(region_from == region_to) ||
            error("Region_From and Region_To cannot be the same. PRAS only considers inter-regional lines in Zonal analysis")

        all([forward_cap, backward_cap] .>= 0.0) ||
            error("Check the forward/backward capacity of Line passed")

        legacy in ["Existing","New"] ||
            error("Unidentified legacy passed")

        0.0 <= FOR <= 1.0 ||
            error("FOR value passed is not allowed")

        MTTR > 0 ||
            error("MTTR value passed is not allowed")

        all([region_from, region_to] .∈  Ref(keys(converter_capacity))) ||
            error("Check the keys of converter capacity dictionary for VSC DC line")

        return new(name, N, category, region_from, region_to, forward_cap, backward_cap, legacy, FOR, MTTR, VSC, converter_capacity)
    end
end

const Lines = Vector{Line}

get_name(ln::Line) = ln.name

get_category(ln::Line) = ln.category

get_forward_capacity(ln::Line) = fill(round(Int,ln.forward_cap),1,ln.N)

get_backward_capacity(ln::Line) = fill(round(Int,ln.backward_cap),1,ln.N)

get_region_from(ln::Line) = ln.region_from

get_region_to(ln::Line) = ln.region_to

get_outage_rate(ln::Line) = outage_to_rate(ln.FOR,ln.MTTR)

get_λ(ln::Line) = fill(getfield(outage_to_rate(ln.FOR,ln.MTTR),:λ),1,ln.N)

get_μ(ln::Line) = fill(getfield(outage_to_rate(ln.FOR,ln.MTTR),:μ),1,ln.N)

function get_legacy_lines(lns::Lines, leg::String)
    if ~(leg in ["Existing","New"])
        error("Unidentified legacy passed")
    end

    leg_line_idxs = findall(getfield.(lns,:legacy) .== leg)
    if isnothing(leg_line_idxs)
        @warn "No lines with this legacy"
    else
        return lns[leg_line_idxs]
    end
end

# Functions for  processing ReEDS2PRAS lines (preparing PRAS lines)
function get_sorted_region_tuples(lns::Lines, region_names::Vector{String})
    regions_from = get_region_from.(lns)
    regions_to = get_region_to.(lns)

    regions_tuple= []
    for (idx,reg) in enumerate(regions_from)
        region_from_idx = findfirst(x->x==reg,region_names)
        region_to_idx = findfirst(x->x==regions_to[idx],region_names)
        if (region_from_idx < region_to_idx)
            push!(regions_tuple,(reg,regions_to[idx]))
        else
            push!(regions_tuple,(regions_to[idx],reg))
        end
    end
    return regions_tuple
end

function get_sorted_region_tuples(lns::Lines)
    regions_from = get_region_from.(lns)
    regions_to = get_region_to.(lns)

    region_names = unique(append!(regions_from, regions_to))

    get_sorted_region_tuples(lns, region_names)
end

function get_sorted_region_tuples(lns::Lines, regions::Vector{Region})
    get_sorted_region_tuples(lns,get_name.(regions))
end

function get_sorted_lines(lns::Lines, region_names::Vector{String})
    regions_tuple = get_sorted_region_tuples(lns, region_names)
    temp_regions_tuple = unique(regions_tuple);
    interface_dict = Dict();

    for i in eachindex(temp_regions_tuple)
        temp = findall(x -> x == temp_regions_tuple[i], regions_tuple);
        push!(interface_dict, temp_regions_tuple[i] => (temp,length(temp)))
    end

    num_interfaces = length(temp_regions_tuple);
    sorted_regional_lines = Line[];
    interface_line_idxs = Array{UnitRange{Int64},1}(undef,num_interfaces);
    start_id = Array{Int64}(undef,num_interfaces); 
    for i in 1: num_interfaces
        for j in interface_dict[temp_regions_tuple[i]][1]
            push!(sorted_regional_lines, lns[j])
        end
        i==1 ? start_id[i] = 1 : start_id[i] =start_id[i-1]+interface_dict[temp_regions_tuple[i-1]][2]
        interface_line_idxs[i] = range(start_id[i], length=interface_dict[temp_regions_tuple[i]][2])
    end

    return sorted_regional_lines, temp_regions_tuple , interface_line_idxs
end

get_sorted_lines(lns::Lines, regions::Vector{Region}) = get_sorted_lines(lns,get_name.(regions))

function process_vsc_lines(lns::Lines, regions::Vector{Region})
    N = first(regions).N
    vsc_dc_line_idxs = findall(getfield.(lns,:VSC))
    vsc_dc_lines_copy = lns[vsc_dc_line_idxs]
    deleteat!(lns,vsc_dc_line_idxs)

    for vsc_line in vsc_dc_lines_copy
        dc_region_from = "DC_$(vsc_line.region_from)"
        dc_region_to = "DC_$(vsc_line.region_to)"

        for reg_name in [dc_region_from, dc_region_to]
            if ~(reg_name in get_name.(regions))
                push!(regions,region(reg_name,N, zeros(Float64,N)))
            end
        end

        push!(lns, line("$(vsc_line.name)_DC", vsc_line.N, vsc_line.category, dc_region_from, dc_region_to, vsc_line.forward_cap,vsc_line.backward_cap,
                        vsc_line.legacy,vsc_line.FOR,vsc_line.MTTR))
        push!(lns, line("$(vsc_line.name)_Converter_From", vsc_line.N, vsc_line.category, dc_region_from, vsc_line.region_from, 
                        vsc_line.converter_capacity[vsc_line.region_from],vsc_line.converter_capacity[vsc_line.region_from],vsc_line.legacy,vsc_line.FOR,vsc_line.MTTR))
        push!(lns, line("$(vsc_line.name)_Converter_To", vsc_line.N, vsc_line.category, dc_region_to, vsc_line.region_to, 
                         vsc_line.converter_capacity[vsc_line.region_to],vsc_line.converter_capacity[vsc_line.region_to],vsc_line.legacy,vsc_line.FOR,vsc_line.MTTR))
        
    end
    return lns,regions
end

function make_pras_interfaces(sorted_lines::Vector{Line},temp_regions_tuples::Vector{Any},interface_line_idxs::Vector{UnitRange{Int64}}, regions::Vector{Region})
    make_pras_interfaces(sorted_lines,temp_regions_tuples,interface_line_idxs, get_name.(regions))
end

function make_pras_interfaces(sorted_lines::Vector{Line},temp_regions_tuples::Vector{Any},interface_line_idxs::Vector{UnitRange{Int64}},region_names::Vector{String})
    num_interfaces = length(temp_regions_tuples);
    interface_regions_from = [findfirst(x->x==temp_regions_tuples[i][1],region_names) for i in 1:num_interfaces];
    interface_regions_to = [findfirst(x->x==temp_regions_tuples[i][2],region_names) for i in 1:num_interfaces];

    N = first(sorted_lines).N
    
    interface_forward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
    interface_backward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);

    line_forward_capacity_array = reduce(vcat,get_forward_capacity.(sorted_lines))
    line_backward_capacity_array = reduce(vcat,get_backward_capacity.(sorted_lines))

    for i in 1:num_interfaces
        interface_forward_capacity_array[i,:] =  sum(line_forward_capacity_array[interface_line_idxs[i],:],dims=1)
        interface_backward_capacity_array[i,:] =  sum(line_backward_capacity_array[interface_line_idxs[i],:],dims=1)
    end

    new_interfaces = PRAS.Interfaces{N,PRAS.MW}(interface_regions_from, interface_regions_to, interface_forward_capacity_array, interface_backward_capacity_array);

    return new_interfaces
end