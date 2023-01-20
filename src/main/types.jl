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
    function Region(name, N, load = zeros(Float64,N))

        0 < N <= 8784 ||
            error("Check the PRAS timesteps (N) passed.")

        length(load) == N ||
            error("The length of the region load time series data should be equal to PRAS timesteps (N)")

        all(load .>= 0.0) ||
            error("Check for inconsistencies in load time series data")

        return new(name, N, load)

    end
end

get_name(reg::Region) = reg.name

get_load(reg::Region) =  permutedims(round.(Int,reg.load))

# Generators
abstract type Generator end
# const Generators = Vector{<:Generator}

struct Thermal_Gen <:Generator
    name::String 
    N::Int64
    region_name::String
    capacity::Float64
    fuel::String
    legacy::String
    FOR::Float64
    MTTR::Int64

    # Inner Constructors & Checks

    function Thermal_Gen(name, N, region_name, capacity, fuel = "OT", legacy = "New", FOR = 0.0, MTTR = 24)

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

    # Inner Constructors & Checks
    function VG_Gen(name, N, region_name, installed_capacity, capacity, type, legacy = "New", FOR = 0.0, MTTR = 24)

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

# Storages
abstract type Storage end

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

    # Inner Constructors & Checks
    function Battery(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, legacy = "New", charge_eff = 1.0, discharge_eff = 1.0, carryover_eff = 1.0,
                     FOR = 0.0, MTTR = 24)

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

    # Inner Constructors & Checks
    function Gen_Storage(name, N, region_name, type, charge_cap, discharge_cap, energy_cap, inflow, grid_withdrawl_cap, grid_inj_cap, legacy = "New", 
        charge_eff = 1.0, discharge_eff = 1.0, carryover_eff = 1.0, FOR = 0.0, MTTR = 24)

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

# Helper Functions
get_outage_rate(stor::STOR) where {STOR <: Storage} = outage_to_rate(stor.FOR,stor.MTTR)

get_λ(stor::STOR) where {STOR <: Storage} = fill(getfield(get_outage_rate(stor),:λ),1,stor.N)

get_μ(stor::STOR) where {STOR <: Storage} = fill(getfield(get_outage_rate(stor),:μ),1,stor.N)

get_category(stor::STOR) where {STOR <: Storage} = "$(stor.legacy)_$(stor.type)"

function get_storages_in_region(stors::Vector{<:Storage}, reg_name::String)
    reg_stors = filter(stor -> stor.region_name == reg_name, stors)
    if isempty(reg_stors)
        # @warn "No storages in region: $(reg_name)"
        return Storage[]
    else
        return reg_stors
    end
end

get_storages_in_region(stors::Vector{<:Storage}, reg::Region) = get_storages_in_region(stors, reg.name)

function get_legacy_storages(stors::Vector{<:Storage}, leg::String)
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

get_components(comps::Vector{<:Generator}, region_name::String) = get_generators_in_region(comps, region_name)
    
emptyvec(::Vector{<:Storage}) = Storage[]

get_components(comps::Vector{<:Storage}, region_name::String) = get_storages_in_region(comps, region_name)

# Functions for processing ReEDS2PRAS generators and storages (preparing PRAS lines)
function get_sorted_components(comps::COMPONENTS, region_names::Vector{String}) where {COMPONENTS <: Union{Vector{<:Generator},Vector{<:Storage}}}
    num_regions = length(region_names)
    all_comps = []
    start_idx = Array{Int64}(undef,num_regions); 
    region_comp_idxs = Array{UnitRange{Int64},1}(undef,num_regions); 

    for (idx,region_name) in enumerate(region_names)
        region_comps = get_components(comps,region_name)
        push!(all_comps,region_comps)
        if idx == 1
            start_idx[idx] = 1
        else
            prev_idx = start_idx[idx-1]
            prev_length = length(all_comps[idx-1])
            start_idx[idx] = prev_idx + prev_length
        end
        region_comp_idxs[idx] = range(start_idx[idx], length=length(all_comps[idx]))
    end
    
    sorted_comps = emptyvec(comps)
    for idx in eachindex(all_comps)
        if (length(all_comps[idx]) != 0)
            append!(sorted_comps,all_comps[idx])
        end
    end
    return sorted_comps, region_comp_idxs
end

get_sorted_components(comps::COMPONENTS, regions::Vector{Region}) where {COMPONENTS <: Union{Vector{<:Generator},Vector{<:Storage}}} = get_sorted_components(comps,get_name.(regions))

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
    
    # Inner Constructors & Checks
    function Line(name, N, category, region_from, region_to, forward_cap, backward_cap, legacy = "New", FOR = 0.0, MTTR = 24, VSC = false, 
        converter_capacity = Dict(region_from => 0.0, region_to => 0.0))

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

get_name(ln::Line) = ln.name

get_category(ln::Line) = ln.category

get_forward_capacity(ln::Line) = fill(round(Int,ln.forward_cap),1,ln.N)

get_backward_capacity(ln::Line) = fill(round(Int,ln.backward_cap),1,ln.N)

get_region_from(ln::Line) = ln.region_from

get_region_to(ln::Line) = ln.region_to

#Helper functions
get_outage_rate(ln::Line) = outage_to_rate(ln.FOR,ln.MTTR)

get_λ(ln::Line) = fill(getfield(outage_to_rate(ln.FOR,ln.MTTR),:λ),1,ln.N)

get_μ(ln::Line) = fill(getfield(outage_to_rate(ln.FOR,ln.MTTR),:μ),1,ln.N)

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
function get_sorted_region_tuples(lines::Vector{Line}, region_names::Vector{String})
    region_idxs = Dict(name => idx for (idx, name) in enumerate(region_names))

    line_from_to_reg_idxs = similar(lines, Tuple{Int,Int})

    for (l, line) in enumerate(lines)

        from_name = get_region_from(line)
        to_name = get_region_to(line)

        from_idx = region_idxs[from_name]       
        to_idx = region_idxs[to_name]

        line_from_to_reg_idxs[l] = from_idx < to_idx ? (from_idx, to_idx) : (to_idx, from_idx)
    end

    return line_from_to_reg_idxs
end

function get_sorted_region_tuples(lines::Vector{Line})
    regions_from = get_region_from.(lines)
    regions_to = get_region_to.(lines)

    region_names = unique(append!(regions_from, regions_to))

    get_sorted_region_tuples(lines, region_names)
end

function get_sorted_region_tuples(lines::Vector{Line}, regions::Vector{Region})
    get_sorted_region_tuples(lines,get_name.(regions))
end

function get_sorted_lines(lines::Vector{Line}, region_names::Vector{String})
    line_from_to_reg_idxs = get_sorted_region_tuples(lines, region_names)
    line_ordering = sortperm(line_from_to_reg_idxs)
    
    sorted_lines = lines[line_ordering]
    sorted_from_to_reg_idxs = line_from_to_reg_idxs[line_ordering] 
    interface_reg_idxs = unique(sorted_from_to_reg_idxs)

    # Ref tells Julia to use interfaces as Vector, only broadcasting over lines_sorted
    interface_line_idxs = searchsorted.(Ref(sorted_from_to_reg_idxs), interface_reg_idxs)

    return sorted_lines, interface_reg_idxs , interface_line_idxs
end

get_sorted_lines(lines::Vector{Line}, regions::Vector{Region}) = get_sorted_lines(lines,get_name.(regions))

function process_vsc_lines(lines::Vector{Line}, regions::Vector{Region})
    N = first(regions).N
    non_vsc_dc_lines = filter(line -> ~line.VSC, lines)
    vsc_dc_lines = filter(line -> line.VSC, lines)
    
    for vsc_line in vsc_dc_lines
        dc_region_from = "DC_$(vsc_line.region_from)"
        dc_region_to = "DC_$(vsc_line.region_to)"

        for reg_name in [dc_region_from, dc_region_to]
            if ~(reg_name in get_name.(regions))
                push!(regions,Region(reg_name,N, zeros(Float64,N)))
            end
        end

        push!(non_vsc_dc_lines, Line("$(vsc_line.name)_DC", vsc_line.N, vsc_line.category, dc_region_from, dc_region_to, vsc_line.forward_cap,vsc_line.backward_cap,
                        vsc_line.legacy,vsc_line.FOR,vsc_line.MTTR))
        push!(non_vsc_dc_lines, Line("$(vsc_line.name)_Converter_From", vsc_line.N, vsc_line.category, dc_region_from, vsc_line.region_from, 
                        vsc_line.converter_capacity[vsc_line.region_from],vsc_line.converter_capacity[vsc_line.region_from],vsc_line.legacy,vsc_line.FOR,vsc_line.MTTR))
        push!(non_vsc_dc_lines, Line("$(vsc_line.name)_Converter_To", vsc_line.N, vsc_line.category, dc_region_to, vsc_line.region_to, 
                         vsc_line.converter_capacity[vsc_line.region_to],vsc_line.converter_capacity[vsc_line.region_to],vsc_line.legacy,vsc_line.FOR,vsc_line.MTTR))
        
    end
    return non_vsc_dc_lines,regions
end

function make_pras_interfaces(sorted_lines::Vector{Line},interface_reg_idxs::Vector{Tuple{Int64, Int64}},interface_line_idxs::Vector{UnitRange{Int64}},
    regions::Vector{Region})
    make_pras_interfaces(sorted_lines,interface_reg_idxs,interface_line_idxs, get_name.(regions))
end

function make_pras_interfaces(sorted_lines::Vector{Line},interface_reg_idxs::Vector{Tuple{Int64, Int64}},interface_line_idxs::Vector{UnitRange{Int64}},
    region_names::Vector{String})

    num_interfaces = length(interface_reg_idxs);
    interface_regions_from = first.(interface_reg_idxs);
    interface_regions_to = last.(interface_reg_idxs);

    N = first(sorted_lines).N

    # Lines
    line_names = get_name.(sorted_lines)
    line_cats = get_category.(sorted_lines)
    line_forward_cap = reduce(vcat,get_forward_capacity.(sorted_lines))
    line_backward_cap = reduce(vcat,get_backward_capacity.(sorted_lines))
    line_λ = reduce(vcat,get_λ.(sorted_lines))
    line_μ = reduce(vcat,get_μ.(sorted_lines))

    new_lines = PRAS.Lines{N,1,PRAS.Hour,PRAS.MW}(line_names, line_cats, line_forward_cap, line_backward_cap, line_λ ,line_μ);

    interface_forward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
    interface_backward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);

    for i in 1:num_interfaces
        interface_forward_capacity_array[i,:] =  sum(line_forward_cap[interface_line_idxs[i],:],dims=1)
        interface_backward_capacity_array[i,:] =  sum(line_backward_cap[interface_line_idxs[i],:],dims=1)
    end

    new_interfaces = PRAS.Interfaces{N,PRAS.MW}(interface_regions_from, interface_regions_to, interface_forward_capacity_array, interface_backward_capacity_array);

    return new_lines, new_interfaces

end