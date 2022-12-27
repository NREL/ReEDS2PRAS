#####################################################
# Surya
# NREL
# December 2022
# ReEDS2PRAS - NTPS
# Structs to make objects from ReEDS data
#######################################################
# Converting FOR and MTTR to λ and μ
#######################################################
function outage_to_rate(outage_data::Tuple{Float64, Int64})
    for_gen = outage_data[1]
    mttr = outage_data[2]

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
struct region
    name::String
    N::Int64
    load::Vector{Int64}

     # Inner Constructors
    region(name = "region_1", N =10) = region(name,N, zeros(Int64,N))

    # Checks
    region(x,y,z) =
    if ~(0 < y <= 8784)
        error("Check the PRAS timesteps (N) passed.")
    elseif (length(z) !== y)
        error("The length of the region load time series data should be equal to PRAS timesteps (N)")
    elseif ~(all(z .>= 0.0))
        "Check for inconsistencies in load time series data"
    else
        new(x,y,z)
    end
end
function get_name(reg::region)
    return reg.name
end

function get_load(reg::region)
    return permutedims(round.(Int,reg.load))
end

# Generators
abstract type generator end
const generators = Vector{<:generator}

struct thermal_gen <:generator
    name::String 
    N::Int64
    region::String
    cap::Float64
    fuel::String
    legacy::String
    FOR::Float64
    MTTR::Int64

    # Inner Constructors
    thermal_gen(s,t,u,v,w,x,y) = thermal_gen(s,t,u,v,w,x,y,24);
    thermal_gen(s,t,u,v,w,x) = thermal_gen(s,t,u,v,w,x,0.0,24);
    thermal_gen(s,t,u,v,w) = thermal_gen(s,t,u,v,w,"New",0.0,24);
    thermal_gen(s,t,u,v) = thermal_gen(s,t,u,v,"OT","New",0.0,24);
    # For illustration purposes
    thermal_gen(nothing) = thermal_gen("thermal_gen_1",10,"reg_1",10.0,"NG","New",0.1,24);

    # Checks
    thermal_gen(s,t,u,v,w,x,y,z) = 
    if ~((0.0 <= y <= 1.0) || (z<=0))
        error("FOR and/or MTTR values passed are not allowed")
    elseif ~(x in ["Existing","New"])
        error("Unidentified legacy passed")
    elseif ~(v >= 0.0)
        error("Generator capacity passed is not allowed") 
    elseif ~(0 < t <= 8784)
        error("Check the PRAS timesteps (N) passed.")
    else
        new(s,t,u,v,w,x,y,z)
    end
end

struct vg_gen <:generator
    name::String 
    N::Int64
    region::String
    installed_cap::Float64
    cap::Vector{Float64}
    type::String
    legacy::String
    FOR::Float64
    MTTR::Int64

    # Inner Constructors
    vg_gen(r,s,t,u,v,w,x) = vg_gen(r,s,t,u,v,w,x,0.0,24)
    vg_gen(r,s,t,u,v,w) =  vg_gen(r,s,t,u,v,w,"New",0.0,24)
    # For illustration purposes
    vg_gen(nothing) = vg_gen("vg_gen_1",10,"reg_1",10.0,zeros(Float64,10),"dupv","New",0.1,24)
   
    # Checks
    vg_gen(r,s,t,u,v,w,x,y,z) =
    if ~((0.0 <= y <= 1.0) || (z<=0))
        error("FOR and/or MTTR values passed are not allowed")
    elseif ~(x in ["Existing","New"])
        error("Unidentified legacy passed")
    elseif ~(0 < s <= 8784)
        error("Check the PRAS timesteps (N) passed.")
    elseif ~(all(0.0 .<= v .<= u))
        error("Check for inconsistencies in VG time series data")
    elseif (length(v) !== s)
        error("The length of the VG time series data should be equal to PRAS timesteps (N)")
    # elseif ~(w in ["wind-ons","wind-ofs","dupv","upv","csp","distpv"])
    #     error("Check the type of VG being passed")
    else 
        new(r,s,t,u,v,w,x,y,z)
    end
end

function get_name(gen::GEN) where {GEN <: generator}
    return gen.name
end

function get_region(gen::GEN) where {GEN <: generator}
    return gen.region
end

function get_nameplate(gen::GEN) where {GEN <: generator}
    return gen.cap
end

function get_capacity(gen::thermal_gen)
    return fill(round(Int,gen.cap),1,gen.N)
end

function get_capacity(gen::vg_gen)
    return permutedims(round.(Int,gen.cap))
end

function get_legacy(gen::GEN) where {GEN <: generator}
    return gen.legacy
end

function get_outage_rate(gen::GEN) where {GEN <: generator}
    rate = outage_to_rate((gen.FOR,gen.MTTR))
    return rate
end

function get_λ(gen::GEN) where {GEN <: generator}
    return fill(getfield(outage_to_rate((gen.FOR,gen.MTTR)),:λ),1,gen.N)
end

function get_μ(gen::GEN) where {GEN <: generator}
    return fill(getfield(outage_to_rate((gen.FOR,gen.MTTR)),:μ),1,gen.N)
end

function get_type(gen::vg_gen)
    return gen.type
end

function get_type(gen::thermal_gen)
    return gen.fuel
end

function get_category(gen::vg_gen)
    return gen.legacy*"_"*gen.type
end

function get_category(gen::thermal_gen)
    return gen.legacy*"_"*"Thermal"*"_"*gen.fuel
end

function get_generators_in_region(gens::generators, reg_name::String)
    reg_gen_idxs = findall(getfield.(gens,:region) .== reg_name)
    if isnothing(reg_gen_idxs)
        @warn "No generators in the region"
    else
        return gens[reg_gen_idxs]
    end
end

function get_generators_in_region(gens::generators, reg::region)
    get_generators_in_region(gens, reg.name)
end

function get_legacy_generators(gens::generators, leg::String)
    if ~(leg in ["Existing","New"])
        error("Unidentified legacy passed")
    end

    leg_gen_idxs = findall(getfield.(gens,:legacy) .== leg)
    if isnothing(leg_gen_idxs)
        @warn "No generators with this legacy"
    else
        return gens[leg_gen_idxs]
    end
end

# Storages
abstract type storage end
const storages = Vector{<:storage}

struct battery <:storage
    name::String 
    N::Int64
    region::String
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
    # battery(name = "stor_1", N = 10, rg = "reg_1",type = "4-hour", c_cap = 10.0, dis_cap = 10.0, energy_cap = 40.0,  lg = "New", chr_eff = 0.9, dis_eff = 1.0, 
    #        cry_eff = 1.0, f_or = 0.0, mttr = 24) = battery(name,N,rg,type, c_cap,dis_cap,energy_cap, lg, chr_eff, dis_eff,cry_eff, f_or, mttr)
    battery(n,o,p,q,r,s,t,u,v,w,x) = battery(n,o,p,q,r,s,t,u,v,w,x,0.0,24)
    battery(n,o,p,q,r,s,t,u,v,w) = battery(n,o,p,q,r,s,t,u,v,w,1.0,0.0,24)
    battery(n,o,p,q,r,s,t,u) = battery(n,o,p,q,r,s,t,u,1.0,1.0,1.0,0.0,24)
    battery(n,o,p,q,r,s,t) = battery(n,o,p,q,r,s,t,"New",1.0,1.0,1.0,0.0,24)
    # For illustration purposes
    battery(nothing) = battery("stor_1",10,"reg_1","4-hour",10.0,10.0,40.0,"New",0.9,1.0,1.0,0.0,24)
   
    # Checks
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) =
    if ~((0.0 <= y <= 1.0) || (z < 0))
        error("FOR and/or MTTR values passed are not allowed")
    elseif ~(u in ["Existing","New"])
        error("Unidentified legacy passed")
    elseif ~(0 < o <= 8784)
        error("Check the PRAS timesteps (N) passed.")
    elseif ~(all(0.0 .<= [v,w,x] .<= 1.0))
        error("Invalid charge/discharge/carryover efficiency passed")
    elseif ~(r > 0.0)
        error("Charge capacity passed is not allowed")
    elseif ~(s > 0.0)
        error("Discharge capacity passed is not allowed")
    elseif ~(t > 0.0)
        error("Energy capacity passed is not allowed")
    else
        new(n,o,p,q,r,s,t,u,v,w,x,y,z)
    end 
end

struct gen_storage <:storage
    name::String 
    N::Int64
    region::String
    type::String
    charge_cap::Vector{Float64}
    discharge_cap::Vector{Float64}
    energy_cap::Vector{Float64}
    inflow::Vector{Float64}
    grid_withdrawl_cap::Vector{Float64}
    grid_inj_cap::Vector{Float64}
    legacy::String
    charge_eff::Float64
    discharge_eff::Float64
    carryover_eff::Float64
    FOR::Float64
    MTTR::Int64

    # Inner Constructors
    # gen_storage(name = "gen_stor_1", N = 10, rg = "reg_1",type = "pumped-storage", c_cap = fill(10.0,N), dis_cap = fill(10.0,N), energy_cap = fill(40.0,N), infl = fill(10.0,N), 
    #             g_with_cap = fill(10.0,N),g_inj_cap = fill(10.0,N), lg = "New", chr_eff = 0.9, dis_eff = 1.0, cry_eff = 1.0, f_or = 0.0, mttr = 24) = 
    #             gen_storage(name,N,rg,type, c_cap,dis_cap,energy_cap, infl, g_with_cap, g_inj_cap, lg, chr_eff, dis_eff,cry_eff, f_or, mttr)

    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x) = gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,0.0,24)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w) = gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,1.0,0.0,24)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u) = gen_storage(k,l,m,n,o,p,q,r,s,t,u,1.0,1.0,1.0,0.0,24)
    gen_storage(k,l,m,n,o,p,q,r,s,t) = gen_storage(k,l,m,n,o,p,q,r,s,t,"New",1.0,1.0,1.0,0.0,24)
    # For illustration purposes
    gen_storage(nothing) = gen_storage("gen_stor_1",10,"reg_1","4-Hour",fill(10.0,10),fill(10.0,10),fill(40.0,10),fill(10.0,10),fill(10.0,10),fill(10.0,10),"New",0.9,1.0,1.0,0.0,24)
   
    # Checks
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = 
    if ~((0.0 <= y <= 1.0) || (z<0))
        error("FOR and/or MTTR values passed are not allowed")
    elseif ~(u in ["Existing","New"])
        error("Unidentified legacy passed")
    elseif ~(0 < l <= 8784)
        error("Check the PRAS timesteps (N) passed.")
    elseif ~(all(0.0 .<= [v,w,x] .<= 1.0))
        error("Invalid charge/discharge/carryover efficiency passed")
    elseif ~(all(length.([o,p,q,r,s,t]) .== l))
        error("The length of the time series data associated with the storage should be equal to PRAS timesteps (N)")
    elseif ~(all(o .>= 0.0))
        error("Check for inconsistencies in generatorstorage charge capacity time series data")
    elseif ~(all(p .>= 0.0))
        error("Check for inconsistencies in generatorstorage discharge capacity time series data")
    elseif ~(all(q .>= 0.0))
        error("Check for inconsistencies in generatorstorage energy capacity time series data")
    elseif ~(all(r .>= 0.0))
        error("Check for inconsistencies in generatorstorage inflow time series data")
    elseif ~(all(s .>= 0.0))
        error("Check for inconsistencies in generatorstorage grid withdrawl capacity time series data")
    elseif ~(all(t .>= 0.0))
        error("Check for inconsistencies in generatorstorage grid injection capacity time series data")
    else
        new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    end
end

function get_name(stor::STOR) where {STOR <: storage}
    return stor.name
end

function get_type(stor::STOR) where {STOR <: storage}
    return stor.type
end

function get_legacy(stor::STOR) where {STOR <: storage}
    return stor.legacy
end

function get_outage_rate(stor::STOR) where {STOR <: storage}
    rate = outage_to_rate((stor.FOR,stor.MTTR))
    return rate
end

function get_λ(stor::STOR) where {STOR <: storage}
    return fill(getfield(outage_to_rate((stor.FOR,stor.MTTR)),:λ),1,stor.N)
end

function get_μ(stor::STOR) where {STOR <: storage}
    return fill(getfield(outage_to_rate((stor.FOR,stor.MTTR)),:μ),1,stor.N)
end

function get_category(stor::STOR) where {STOR <: storage}
    return stor.legacy*"_"*stor.type
end

function get_charge_capacity(stor::battery)
    return fill(round(Int,stor.charge_cap),1,stor.N)
end

function get_charge_capacity(stor::gen_storage)
    return permutedims(round.(Int,stor.charge_cap))
end

function get_discharge_capacity(stor::battery)
    return fill(round(Int,stor.discharge_cap),1,stor.N)
end

function get_discharge_capacity(stor::gen_storage)
    return permutedims(round.(Int,stor.discharge_cap))
end

function get_energy_capacity(stor::battery)
    return fill(round(Int,stor.energy_cap),1,stor.N)
end

function get_energy_capacity(stor::gen_storage)
    return permutedims(round.(Int,stor.energy_cap))
end

function get_inflow(stor::gen_storage)
    return permutedims(round.(Int,stor.inflow))
end

function get_grid_withdrawl_capacity(stor::gen_storage)
    return permutedims(round.(Int,stor.grid_withdrawl_cap))
end

function get_grid_injection_capacity(stor::gen_storage)
    return permutedims(round.(Int,stor.grid_inj_cap))
end

function get_charge_efficiency(stor::STOR) where {STOR <: storage}
    return fill(stor.charge_eff,1,stor.N)
end

function get_discharge_efficiency(stor::STOR) where {STOR <: storage}
    return fill(stor.discharge_eff,1,stor.N)
end

function get_carryover_efficiency(stor::STOR) where {STOR <: storage}
    return fill(stor.carryover_eff,1,stor.N)
end

function get_storages_in_region(stors::storages, reg_name::String)
    reg_stor_idxs = findall(getfield.(stors,:region) .== reg_name)
    if isnothing(reg_stor_idxs)
        @warn "No storages in the region"
    else
        return stors[reg_stor_idxs]
    end
end

function get_storages_in_region(stors::storages, reg::region)
    get_storages_in_region(stors, reg.name)
end

function get_legacy_storages(stors::storages, leg::String)
    if ~(leg in ["Existing","New"])
        error("Unidentified legacy passed")
    end

    leg_stor_idxs = findall(getfield.(stors,:legacy) .== leg)
    if isnothing(eg_stor_idxs)
        @warn "No storages with this legacy"
    else
        return stors[leg_stor_idxs]
    end
end
component_dict = Dict(
    Vector{generator} => (func = get_generators_in_region, vec = generator[]),
    Vector{storage} => (func = get_storages_in_region, vec = storage[]),
    Vector{gen_storage} => (func = get_storages_in_region, vec = storage[])
)
# Functions for processing ReEDS2PRAS generators and storages (preparing PRAS lines)
function get_sorted_components(comps::COMPONENTS, region_names::Vector{String}) where {COMPONENTS <: Union{generators,storages}}
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

function get_sorted_components(comps::COMPONENTS, regions::Vector{region}) where {COMPONENTS <: Union{generators,storages}}
    get_sorted_components(comps,get_name.(regions))
end

# Lines
struct line
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

    # Inner Constructors
    # line(name = "line_1", N = 10, cat = "AC", reg_from = "1", reg_to = "2", for_cap = 10.0, back_cap = 10.0, leg = "Existing",f_or = 0.0, mttr = 24) = 
    #      line(name, N, cat,reg_from,reg_to, for_cap, back_cap, leg, f_or, mttr)
    line(q,r,s,t,u,v,w,x) = line(q,r,s,t,u,v,w,x,0.0,24)
    line(q,r,s,t,u,v,w) = line(q,r,s,t,u,v,w,"New",0.0,24)
    line(q,r,s,t,u,v) = line(q,r,s,t,u,v,v,"New",0.0,24)
    # For illustration purposes
    line(nothing) = line("line_1",10,"AC","1","2",10.0,10.0,"New",0.0,24)

    # Checks
    line(q,r,s,t,u,v,w,x,y,z) = 
    if ~(0 < r <= 8784)
        error("Check the PRAS timesteps (N) passed.")
    elseif ~(s in ["AC","Two-Terminal","VSC-DC"])
        error("Check the category of line passed")
    elseif (t == u)
        error("Region From and Region To cannot be the same. PRAS only considers inter-regional lines in Zonal analysis")
    elseif ~(all([v,w] .> 0.0))
        error("Check the forward/backward capacity of line passed")
    elseif ~(x in ["Existing","New"])
        error("Unidentified legacy passed")
    elseif ~((0.0 <= y <= 1.0) || (z<0))
        error("FOR and/or MTTR values passed are not allowed")
    else
        new(q,r,s,t,u,v,w,x,y,z)
    end
end

const lines = Vector{line}

function get_name(ln::line)
    return ln.name
end

function get_category(ln::line)
    return ln.category
end

function get_forward_capacity(ln::line)
    return fill(round(Int,ln.forward_cap),1,ln.N)
end

function get_backward_capacity(ln::line)
    return fill(round(Int,ln.backward_cap),1,ln.N)
end

function get_region_from(ln::line)
    return ln.region_from
end

function get_region_to(ln::line)
    return ln.region_to
end

function get_outage_rate(ln::line)
    rate = outage_to_rate((ln.FOR,ln.MTTR))
    return rate
end

function get_λ(ln::line)
    return fill(getfield(outage_to_rate((ln.FOR,ln.MTTR)),:λ),1,ln.N)
end

function get_μ(ln::line)
    return fill(getfield(outage_to_rate((ln.FOR,ln.MTTR)),:μ),1,ln.N)
end

function get_legacy_lines(lns::lines, leg::String)
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
function get_sorted_region_tuples(lns::lines, region_names::Vector{String})
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

function get_sorted_region_tuples(lns::lines)
    regions_from = get_region_from.(lns)
    regions_to = get_region_to.(lns)

    region_names = unique(append!(regions_from, regions_to))

    get_sorted_region_tuples(lns, region_names)
end

function get_sorted_region_tuples(lns::lines, regions::Vector{region})
    get_sorted_region_tuples(lns,get_name.(regions))
end

function get_sorted_lines(lns::lines, region_names::Vector{String})
    regions_tuple = get_sorted_region_tuples(lns, region_names)
    temp_regions_tuple = unique(regions_tuple);
    interface_dict = Dict();

    for i in eachindex(temp_regions_tuple)
        temp = findall(x -> x == temp_regions_tuple[i], regions_tuple);
        push!(interface_dict, temp_regions_tuple[i] => (temp,length(temp)))
    end

    num_interfaces = length(temp_regions_tuple);
    sorted_regional_lines = line[];
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

function get_sorted_lines(lns::lines, regions::Vector{region})
    get_sorted_lines(lns,get_name.(regions))
end

function make_pras_interfaces(sorted_lines::Vector{line},temp_regions_tuples::Vector{Any},interface_line_idxs::Vector{UnitRange{Int64}}, regions::Vector{region})
    make_pras_interfaces(sorted_lines,temp_regions_tuples,interface_line_idxs, get_name.(regions))
end

function make_pras_interfaces(sorted_lines::Vector{line},temp_regions_tuples::Vector{Any},interface_line_idxs::Vector{UnitRange{Int64}},region_names::Vector{String})
    num_interfaces = length(temp_regions_tuple);
    interface_regions_from = [findfirst(x->x==temp_regions_tuples[i][1],region_names) for i in 1:num_interfaces];
    interface_regions_to = [findfirst(x->x==temp_regions_tuples[i][2],region_names) for i in 1:num_interfaces];

    N = first(sorted_lines).N
    
    interface_forward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
    interface_backward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);

    line_forward_capacity_array = reduce(vcat,get_forward_capacity.(sorted_regional_lines))
    line_backward_capacity_array = reduce(vcat,get_backward_capacity.(sorted_regional_lines))

    for i in 1:num_interfaces
        interface_forward_capacity_array[i,:] =  sum(line_forward_capacity_array[interface_line_idxs[i],:],dims=1)
        interface_backward_capacity_array[i,:] =  sum(line_backward_capacity_array[interface_line_idxs[i],:],dims=1)
    end

    new_interfaces = PRAS.Interfaces{N,PRAS.MW}(interface_regions_from, interface_regions_to, interface_forward_capacity_array, interface_backward_capacity_array);

    return new_interfaces
end
# Testing
import PRAS
N=10
# Regions
regs = region[]
push!(regs, region("reg_1",10, fill(100.0,10)))
push!(regs, region("reg_2",10, fill(100.0,10)))

reg_names = get_name.(regs)
reg_load = reduce(vcat,get_load.(regs))

new_regions = PRAS.Regions{N,PRAS.MW}(reg_names, reg_load);
# Generators
gens = generator[]
push!(gens,thermal_gen("thermal_gen_1", 10, "reg_1", 10.0, "NG", "New", 0.1, 24))
push!(gens,thermal_gen("thermal_gen_2", 10, "reg_2", 10.0, "NG", "New", 0.1, 24))
push!(gens,vg_gen("vg_gen_1", 10, "reg_1", 10.0, ones(Float64,10), "dupv", "New", 0.1, 24))
push!(gens,vg_gen("vg_gen_2", 10, "reg_2", 10.0, ones(Float64,10), "dupv", "New", 0.1, 24))

sorted_gens, region_gen_idxs = get_sorted_components(gens,regs)

gen_names = get_name.(sorted_gens)
gen_cats = get_category.(sorted_gens)
gen_cap = reduce(vcat,get_capacity.(sorted_gens))
gen_λ = reduce(vcat,get_λ.(sorted_gens))
gen_μ = reduce(vcat,get_μ.(sorted_gens))

new_generators = PRAS.Generators{N,1,PRAS.Hour,PRAS.MW}(gen_names, gen_cats, gen_cap , gen_λ ,gen_μ);

stors = storage[]
push!(stors,battery("stor_1", 10, "reg_1", "4-hour", 10.0, 10.0, 40.0, "New", 0.9, 1.0, 1.0, 0.0, 24))
push!(stors,battery("stor_2", 10, "reg_2", "4-hour", 10.0, 10.0, 40.0, "New", 0.9, 1.0, 1.0, 0.0, 24))

sorted_stors, reg_stor_idxs  = get_sorted_components(stors,regs)

stor_names = get_name.(sorted_stors)
stor_cats = get_category.(sorted_stors)
stor_cap_array = reduce(vcat,get_charge_capacity.(sorted_stors))
stor_dis_cap_array = reduce(vcat,get_discharge_capacity.(sorted_stors))
stor_enrgy_cap_array = reduce(vcat,get_energy_capacity.(sorted_stors))
stor_chrg_eff_array = reduce(vcat,get_charge_efficiency.(sorted_stors))
stor_dischrg_eff_array = reduce(vcat,get_discharge_efficiency.(sorted_stors))
stor_carryovr_eff_array = reduce(vcat,get_carryover_efficiency.(sorted_stors))
stor_λ = reduce(vcat,get_λ.(sorted_stors))
stor_μ = reduce(vcat,get_μ.(sorted_stors))

new_storage = PRAS.Storages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(stor_names,stor_cats,stor_cap_array,stor_dis_cap_array,stor_enrgy_cap_array,
              stor_chrg_eff_array,stor_dischrg_eff_array, stor_carryovr_eff_array,stor_λ,stor_μ);

# GeneratorStorages
gen_stors = gen_storage[]
push!(gen_stors,gen_storage("gen_stor_1", 10, "reg_1", "Pumped-Hydro", fill(10.0,10),fill(10.0,10), fill(40.0,10),fill(10.0,10),fill(10.0,10),fill(10.0,10),
                            "New", 0.9, 1.0, 1.0, 0.0, 24))
push!(gen_stors,gen_storage("gen_stor_2", 10, "reg_2", "Pumped-Hydro", fill(10.0,10),fill(10.0,10), fill(40.0,10),fill(10.0,10),fill(10.0,10),fill(10.0,10),
                            "New", 0.9, 1.0, 1.0, 0.0, 24))

sorted_gen_stors, reg_genstor_idxs  = get_sorted_components(gen_stors,regs);

gen_stor_names = get_name.(sorted_gen_stors)
gen_stor_cats = get_category.(sorted_gen_stors)
gen_stor_cap_array = reduce(vcat,get_charge_capacity.(sorted_gen_stors))
gen_stor_dis_cap_array = reduce(vcat,get_discharge_capacity.(sorted_gen_stors))
gen_stor_enrgy_cap_array = reduce(vcat,get_energy_capacity.(sorted_gen_stors))
gen_stor_chrg_eff_array = reduce(vcat,get_charge_efficiency.(sorted_gen_stors))
gen_stor_dischrg_eff_array = reduce(vcat,get_discharge_efficiency.(sorted_gen_stors))
gen_stor_carryovr_eff_array = reduce(vcat,get_carryover_efficiency.(sorted_gen_stors))
gen_stor_inflow_array = reduce(vcat,get_inflow.(sorted_gen_stors))
gen_stor_grid_withdrawl_array = reduce(vcat,get_grid_withdrawl_capacity.(sorted_gen_stors))
gen_stor_grid_inj_array = reduce(vcat,get_grid_injection_capacity.(sorted_gen_stors))
gen_stor_λ = reduce(vcat,get_λ.(sorted_gen_stors))
gen_stor_μ = reduce(vcat,get_μ.(sorted_gen_stors))

new_gen_stors = PRAS.GeneratorStorages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(gen_stor_names,gen_stor_cats,gen_stor_cap_array, gen_stor_dis_cap_array, gen_stor_enrgy_cap_array,
                                                                       gen_stor_chrg_eff_array, gen_stor_dischrg_eff_array, gen_stor_carryovr_eff_array,gen_stor_inflow_array,
                                                                       gen_stor_grid_withdrawl_array, gen_stor_grid_inj_array,gen_stor_λ,gen_stor_μ);
all_lines = line[]
push!(all_lines, line("line_1_2", 10, "AC", "reg_1", "reg_2", 10.0, 10.0, "New", 0.0, 24))
push!(all_lines,line("line_2_1", 10, "AC", "reg_2", "reg_1", 10.0, 10.0, "New", 0.0, 24))

sorted_regional_lines,temp_regions_tuple,interface_line_idxs = get_sorted_lines(all_lines,regs)

line_names = get_name.(sorted_regional_lines)
line_cats = get_category.(sorted_regional_lines)
line_forward_cap = reduce(vcat,get_forward_capacity.(sorted_regional_lines))
line_backward_cap = reduce(vcat,get_backward_capacity.(sorted_regional_lines))
line_λ = reduce(vcat,get_λ.(sorted_regional_lines))
line_μ = reduce(vcat,get_μ.(sorted_regional_lines))

new_lines = PRAS.Lines{N,1,PRAS.Hour,PRAS.MW}(line_names, line_cats, line_forward_cap, line_backward_cap, line_λ ,line_μ);
new_interfaces = make_pras_interfaces(sorted_regional_lines,temp_regions_tuple,interface_line_idxs,regs);

import Dates
import TimeZones
first_ts = TimeZones.ZonedDateTime(2023, 01, 01, 00, TimeZones.tz"UTC")
last_ts = first_ts + Dates.Hour(9)
my_timestamps = StepRange(first_ts, Dates.Hour(1), last_ts);

pras_system = PRAS.SystemModel(new_regions, new_interfaces, new_generators, region_gen_idxs, new_storage, reg_stor_idxs, new_gen_stors,
                               reg_genstor_idxs, new_lines,interface_line_idxs,my_timestamps);