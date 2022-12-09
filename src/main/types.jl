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
# TODO : Change load to Int64
struct region
    name::String
    N::Int64
    load::Vector{Float64}

     # Inner Constructors
    region(name = "region_1", N = 10, load = ones(Float64,N)) = new(name,N,load)
    region(name = "region_1", N =10) = new(name,N, zeros(Float64,N))

    # Checks
    region(x,y,z) = !(0 < y <= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(x,y,z)
    region(x,y,z) = (length(z) !== y) ? error("The length of the region load time series data should be equal to PRAS timesteps (N)") : new(x,y,z)
    region(x,y,z) = !(all(z .>= 0.0)) ? error("Check for inconsistencies in load time series data") : new(x,y,z)
end

function get_name(reg::region)
    return reg.name
end

function get_load(reg::region)
    return reg.load
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
    thermal_gen(name = "gen_1", N = 10, rg = "1", cap = 10.0, fl = "NG", lg = "New", f_or = 0.1, mttr = 24) = new(name,N,rg,cap,fl,lg,f_or,mttr)
    thermal_gen(s,t,u,v,w,x,y) = new(s,t,u,v,w,x,y,24)
    thermal_gen(s,t,u,v,w,x) = new(s,t,u,v,w,x,0.0,24)
    thermal_gen(s,t,u,v,w) = new(s,t,u,v,w,"New",0.0,24)
    thermal_gen(s,t,u,v) = new(s,t,u,v,"OT","New",0.0,24)
    # For illustration purposes
    thermal_gen(nothing) = new("thermal_gen_1",10,"reg_1",10.0,"NG","New",0.1,24)

    # Checks
    thermal_gen(s,t,u,v,w,x,y,z) = !((0.0 <= y <= 1.0) || (z<=0)) ? error("FOR and/or MTTR values passed are not allowed") : new(s,t,u,v,w,x,y,z)
    thermal_gen(s,t,u,v,w,x,y,z) = !(x in ["Existing","New"]) ? error("Unidentified legacy passed") : new(s,t,u,v,w,x,y,z)
    thermal_gen(s,t,u,v,w,x,y,z) = !(v > 0.0) ? error("Generator cap passed is not allowed") : new(s,t,u,v,w,x,y,z)
    thermal_gen(s,t,u,v,w,x,y,z) = !(0 < t <= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(s,t,u,v,w,x,y,z)
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
    vg_gen(name = "gen_1", N = 10, rg = "1", inst_cap = 10.0,cap = zeros(Float64,N), type = "dupv", lg = "New", f_or = 0.0, mttr = 24) = 
           new(name,N,rg,inst_cap,cap,type,lg,f_or,mttr)
    vg_gen(r,s,t,u,v,w,x) = new(r,s,t,u,v,w,x,0.0,24)
    vg_gen(r,s,t,u,v,w) = new(r,s,t,u,v,w,"New",0.0,24)
    # For illustration purposes
    vg_gen(nothing) = new("vg_gen_1",10,"reg_1",10.0,zeros(Float64,10),"PV","New",0.1,24)
   
    # Checks
    vg_gen(r,s,t,u,v,w,x,y,z) = !((0.0 <= y <= 1.0) || (z<=0)) ? error("FOR and/or MTTR values passed are not allowed") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(x in ["Existing","New"]) ? error("Unidentified legacy passed") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(0 < s <= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(all(0.0 .<= v .<= u)) ? error("Check for inconsistencies in VG time series data") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(length(v) !== s) ? error("The length of the VG time series data should be equal to PRAS timesteps (N)") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(w in ["wind-ons","wind-ofs","dupv","upv","csp","distpv"]) ? error("Check the type of VG being passed") : new(r,s,t,u,v,w,x,y,z)
end

function get_name(gen::GEN) where {GEN <: generator}
    return gen.name
end

function get_capacity(gen::thermal_gen)
    return fill(round(Int,gen.cap),gen.N)
end

function get_capacity(gen::vg_gen)
    return round.(Int,gen.cap)
end

function get_legacy(gen::GEN) where {GEN <: generator}
    return gen.legacy
end

function get_outage_rate(gen::GEN) where {GEN <: generator}
    rate = outage_to_rate((gen.FOR,gen.MTTR))
    return rate
end

function get_λ(gen::GEN) where {GEN <: generator}
    return fill(getfield(outage_to_rate((gen.FOR,gen.MTTR)),:λ),gen.N)
end

function get_μ(gen::GEN) where {GEN <: generator}
    return fill(getfield(outage_to_rate((gen.FOR,gen.MTTR)),:μ),gen.N)
end

function get_fuel(gen::thermal_gen)
    return gen.fuel
end

function get_type(gen::vg_gen)
    return gen.type
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
    battery(name = "stor_1", N = 10, rg = "reg_1",type = "4-hour", c_cap = 10.0, dis_cap = 10.0, energy_cap = 40.0,  lg = "New", chr_eff = 0.9, dis_eff = 1.0, 
           cry_eff = 1.0, f_or = 0.0, mttr = 24) = new(name,N,rg,type, c_cap,dis_cap,energy_cap, lg, chr_eff, dis_eff,cry_eff, f_or, mttr)
    battery(n,o,p,q,r,s,t,u,v,w,x) = new(n,o,p,q,r,s,t,u,v,w,x,0.0,24)
    battery(n,o,p,q,r,s,t,u,v,w) = new(n,o,p,q,r,s,t,u,v,w,1.0,0.0,24)
    battery(n,o,p,q,r,s,t,u) = new(n,o,p,q,r,s,t,u,1.0,1.0,1.0,0.0,24)
    battery(n,o,p,q,r,s,t) = new(n,o,p,q,r,s,t,"New",1.0,1.0,1.0,0.0,24)
    # For illustration purposes
    battery(nothing) = new("stor_1",10,"reg_1","4-hour",10.0,10.0,40.0,"New",0.9,1.0,1.0,0.0,24)
   
    # Checks
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = (!(0.0 <= y <= 1.0) || (z < 0)) ? error("FOR and/or MTTR values passed are not allowed") : new(n,o,p,q,r,s,t,u,v,w,x,y,z)
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = !(u in ["Existing","New"]) ? error("Unidentified legacy passed") : new(n,o,p,q,r,s,t,u,v,w,x,y,z)
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = !(0 < o <= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(n,o,p,q,r,s,t,u,v,w,x,y,z)
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = !(all(0.0 .<= [v,w,x] .<= 1.0)) ? error("Invalid charge/discharge/carryover efficiency passed") : new(n,o,p,q,r,s,t,u,v,w,x,y,z)
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = !(r > 0.0) ? error("Charge cap passed is not allowed") : new(n,o,p,q,r,s,t,u,v,w,x,y,z) 
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = !(s > 0.0) ? error("Discharge cap passed is not allowed") : new(n,o,p,q,r,s,t,u,v,w,x,y,z) 
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = !(t > 0.0) ? error("Energy capacity passed is not allowed") : new(n,o,p,q,r,s,t,u,v,w,x,y,z) 
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
    gen_storage(name = "gen_stor_1", N = 10, rg = "reg_1",type = "pumped-storage", c_cap = fill(10.0,N), dis_cap = fill(10.0,N), energy_cap = fill(40.0,N), infl = fill(10.0,N), 
                g_with_cap = fill(10.0,N),g_inj_cap = fill(10.0,N), lg = "New", chr_eff = 0.9, dis_eff = 1.0, cry_eff = 1.0, f_or = 0.0, mttr = 24) = 
                new(name,N,rg,type, c_cap,dis_cap,energy_cap, infl, g_with_cap, g_inj_cap, lg, chr_eff, dis_eff,cry_eff, f_or, mttr)

    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x) = new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,0.0,24)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w) = new(k,l,m,n,o,p,q,r,s,t,u,v,w,1.0,0.0,24)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u) = new(k,l,m,n,o,p,q,r,s,t,u,1.0,1.0,1.0,0.0,24)
    gen_storage(k,l,m,n,o,p,q,r,s,t) = new(k,l,m,n,o,p,q,r,s,t,"New",1.0,1.0,1.0,0.0,24)
    # For illustration purposes
    gen_storage(nothing) = new("gen_stor_1",10,"reg_1","pumped-storage",fill(10.0,10),fill(10.0,10),fill(40.0,10),fill(10.0,10),fill(10.0,10),fill(10.0,10),"New",0.9,1.0,1.0,0.0,24)
   
    # Checks
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !((0.0 <= y <= 1.0) || (z<0)) ? error("FOR and/or MTTR values passed are not allowed") : new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(u in ["Existing","New"]) ? error("Unidentified legacy passed") : new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(0 < l <= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(all(0.0 .<= [v,w,x] .<= 1.0)) ? error("Invalid charge/discharge/carryover efficiency passed") : 
                                                   new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(all(length.(o,p,q,r,s,t) .== l)) ? 
                                                   error("The length of the time series data associated with the storage should be equal to PRAS timesteps (N)") :
                                                   new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(all(o .>= 0.0)) ? error("Check for inconsistencies in generatorstorage charge capacity time series data") : 
                                                   new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(all(p .>= 0.0)) ? error("Check for inconsistencies in generatorstorage discharge capacity time series data") : 
                                                   new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(all(q .>= 0.0)) ? error("Check for inconsistencies in generatorstorage energy capacity time series data") : 
                                                   new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(all(r .>= 0.0)) ? error("Check for inconsistencies in generatorstorage inflow time series data") : 
                                                   new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(all(s .>= 0.0)) ? error("Check for inconsistencies in generatorstorage grid withdrawl capacity time series data") : 
                                                   new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(all(t .>= 0.0)) ? error("Check for inconsistencies in generatorstorage grid injection capacity time series data") : 
                                                   new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    
end

function get_name(stor::STOR) where {STOR <: storage}
    return stor.name
end

function get_legacy(stor::STOR) where {STOR <: storage}
    return stor.legacy
end

function get_outage_rate(stor::STOR) where {STOR <: storage}
    rate = outage_to_rate((stor.FOR,stor.MTTR))
    return rate
end

function get_λ(stor::STOR) where {STOR <: storage}
    return fill(getfield(outage_to_rate((stor.FOR,stor.MTTR)),:λ),stor.N)
end

function get_μ(stor::STOR) where {STOR <: storage}
    return fill(getfield(outage_to_rate((stor.FOR,stor.MTTR)),:μ),stor.N)
end

function get_category(stor::STOR) where {STOR <: storage}
    return stor.legacy*"_"*stor.type
end

function get_charge_capacity(stor::battery)
    return fill(round(Int,stor.charge_cap),stor.N)
end

function get_charge_capacity(stor::gen_storage)
    return round.(Int,stor.charge_cap)
end

function get_discharge_capacity(stor::battery)
    return fill(round(Int,stor.discharge_cap),stor.N)
end

function get_discharge_capacity(stor::gen_storage)
    return round.(Int,stor.discharge_cap)
end

function get_energy_capacity(stor::battery)
    return fill(round(Int,stor.energy_cap),stor.N)
end

function get_energy_capacity(stor::gen_storage)
    return round.(Int,stor.energy_cap)
end

function get_inflow(stor::gen_storage)
    return round.(Int,stor.inflow)
end

function get_grid_withdrawl_capacity(stor::gen_storage)
    return round.(Int,stor.grid_withdrawl_cap)
end

function get_grid_injection_capacity(stor::gen_storage)
    return round.(Int,stor.grid_inj_cap)
end

function get_charge_efficiency(stor::STOR) where {STOR <: storage}
    return fill(stor.charge_eff,stor.N)
end

function get_discharge_efficiency(stor::STOR) where {STOR <: storage}
    return fill(stor.discharge_eff,stor.N)
end

function get_carryover_efficiency(stor::STOR) where {STOR <: storage}
    return fill(stor.carryover_eff,stor.N)
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
    leg_stor_idxs = findall(getfield.(stors,:legacy) .== leg)
    if isnothing(eg_stor_idxs)
        @warn "No storages with this legacy"
    else
        return stors[leg_stor_idxs]
    end
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
    line(name = "line_1", N = 10, cat = "AC", reg_from = "1", reg_to = "2", for_cap = 10.0, back_cap = 10.0, leg = "existing",f_or = 0.0, mttr = 24) = 
        new(name, N, cat,reg_from,reg_to, for_cap, back_cap, leg, f_or, mttr)
    line(q,r,s,t,u,v,w,x) = new(q,r,s,t,u,v,w,x,0.0,24)
    line(q,r,s,t,u,v,w) = new(q,r,s,t,u,v,w,"new",0.0,24)
    line(q,r,s,t,u,v) = new(q,r,s,t,u,v,v,"new",0.0,24)
    # For illustration purposes
    line(nothing) = new("line_1",10,"AC","1","2",10.0,10.0,"new",0.0,24)

    # Checks
    line(q,r,s,t,u,v,w,x,y,z)= !(0 < r <= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(q,r,s,t,u,v,w,x,y,z)
    line(q,r,s,t,u,v,w,x,y,z)= !(s in ["AC","Two-Terminal","VSC-DC"]) ? error("Check the category of line passed") : new(q,r,s,t,u,v,w,x,y,z)
    line(q,r,s,t,u,v,w,x,y,z)= (t == u) ? error("region from and region to cannot be the same. PRAS only considers inter-regional lines in Zonal analysis") : 
                               new(q,r,s,t,u,v,w,x,y,z)
    line(q,r,s,t,u,v,w,x,y,z)= !(all([v,w] .> 0.0)) ? error("Check the forward/backward capacity of line passed") : new(q,r,s,t,u,v,w,x,y,z)
    line(q,r,s,t,u,v,w,x,y,z)= !(x in ["Existing","New"]) ? error("Unidentified legacy passed") : new(q,r,s,t,u,v,w,x,y,z)
    line(q,r,s,t,u,v,w,x,y,z)= !((0.0 <= y <= 1.0) || (z<0)) ? error("FOR and/or MTTR values passed are not allowed") : new(q,r,s,t,u,v,w,x,y,z)
end

# TODO implement methods to access fields of line objects

# Testing
gens = generator[]
push!(gens,thermal_gen(nothing))
push!(gens,vg_gen(nothing))

gen_names = get_name.(gens)
gen_cats = get_category.(gens)
gen_cap = get_capacity.(gens)
gen_λ = get_λ.(gens)
gen_μ = get_μ.(gens)

stors = storage[]
push!(stors,battery(nothing))
push!(stors,battery(nothing))

stor_names = get_name.(stors)
stor_cats = get_category.(stors)
stor_cap_array = get_charge_capacity.(stors)
stor_dis_cap_array = get_discharge_capacity.(stors)
stor_enrgy_cap_array = get_energy_capacity.(stors)
stor_chrg_eff_array = get_charge_efficiency.(stors)
stor_dischrg_eff_array = get_discharge_efficiency.(stors)
stor_carryovr_eff_array = get_carryover_efficiency.(stors)
stor_λ = get_λ.(stors)
stor_μ = get_μ.(stors)

gen_stors = storage[]
push!(gen_stors,gen_storage(nothing))
push!(gen_stors,gen_storage(nothing))

gen_stor_names = get_name.(gen_stors)
gen_stor_cats = get_category.(gen_stors)
gen_stor_cap_array = get_charge_capacity.(gen_stors)
gen_stor_dis_cap_array = get_discharge_capacity.(gen_stors)
gen_stor_enrgy_cap_array = get_energy_capacity.(gen_stors)
gen_stor_chrg_eff_array = get_charge_efficiency.(gen_stors)
gen_stor_dischrg_eff_array = get_discharge_efficiency.(gen_stors)
gen_stor_carryovr_eff_array = get_carryover_efficiency.(gen_stors)
gen_stor_inflow_array = get_inflow.(gen_stors)
gen_stor_grid_withdrawl_array = get_grid_withdrawl_capacity.(gen_stors)
gen_stor_grid_inj_array = get_grid_injection_capacity.(gen_stors)
gen_stor_λ = get_λ.(gen_stors)
gen_stor_μ = get_μ.(gen_stors)

lines = line[]
push!(lines, line(nothing))
push!(lines, line(nothing))