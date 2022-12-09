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
    region(x,y,z) = !(0 <= y <= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(x,y,z)
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
    thermal_gen(name = "gen_1", N = 10, rg = "1", cap = 10.0, fl = "NG", lg = "new", f_or = 0.1, mttr = 24) = new(name,N,rg,cap,fl,lg,f_or,mttr)
    thermal_gen(s,t,u,v,w,x,y) = new(s,t,u,v,w,x,y,24)
    thermal_gen(s,t,u,v,w,x) = new(s,t,u,v,w,x,0.0,24)
    thermal_gen(s,t,u,v,w) = new(s,t,u,v,w,"new",0.0,24)
    thermal_gen(s,t,u,v) = new(s,t,u,v,"OT","new",0.0,24)

    # Checks
    thermal_gen(s,t,u,v,w,x,y,z) = !((0.0 <= y <= 1.0) || (z<=0)) ? error("FOR and/or MTTR values passed are not allowed") : new(s,t,u,v,w,x,y,z)
    thermal_gen(s,t,u,v,w,x,y,z) = !(x in ["old","new"]) ? error("Unidentified legacy passed") : new(s,t,u,v,w,x,y,z)
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
    vg_gen(name = "gen_1", N = 10, rg = "1", inst_cap = 10.0,cap = zeros(Float64,N), type = "dupv", lg = "new", f_or = 0.0, mttr = 24) = 
           new(name,N,rg,inst_cap,cap,type,lg,f_or,mttr)
    vg_gen(r,s,t,u,v,w,x) = new(r,s,t,u,v,w,x,0.0,24)
    vg_gen(r,s,t,u,v,w) = new(r,s,t,u,v,w,"new",0.0,24)
   
    # Checks
    vg_gen(r,s,t,u,v,w,x,y,z) = !((0.0 <= y <= 1.0) || (z<=0)) ? error("FOR and/or MTTR values passed are not allowed") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(x in ["old","new"]) ? error("Unidentified legacy passed") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(0< s <= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(all(0.0 .<= v .<= u)) ? error("Check for inconsistencies in VG time series data") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(length(v) !== s) ? error("The length of the VG time series data should be equal to PRAS timesteps (N)") : new(r,s,t,u,v,w,x,y,z)
    vg_gen(r,s,t,u,v,w,x,y,z) = !(w in ["wind-ons","wind-ofs","dupv","upv","csp","distpv"]) ? error("Check the type of VG being passed") : new(r,s,t,u,v,w,x,y,z)
end

function get_name(gen::GEN) where {GEN <: generator}
    return gen.name
end

function get_capacity(gen::thermal_gen)
    return fill(round(Int,gen.cap),1,gen.N)
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
    return getfield(outage_to_rate((gen.FOR,gen.MTTR)),:λ)
end

function get_μ(gen::GEN) where {GEN <: generator}
    return getfield(outage_to_rate((gen.FOR,gen.MTTR)),:μ)
end

function get_fuel(gen::thermal_gen)
    return gen.fuel
end

function get_type(gen::vg_gen)
    return gen.type
end

function get_category(gen::vg_gen)
    return gen.legacy*gen.type
end

function get_category(gen::thermal_gen)
    return gen.legacy*"Thermal"*gen.fuel
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
    battery(name = "stor_1", N = 10, rg = "reg_1",type = "4-hour", c_cap = 10.0, dis_cap = 10.0, energy_cap = 40.0,  lg = "new", chr_eff = 0.9, dis_eff = 1.0, 
           cry_eff = 1.0, f_or = 0.0, mttr = 24) = new(name,N,rg,type, c_cap,dis_cap,energy_cap, lg, chr_eff, dis_eff,cry_eff, f_or, mttr)
    battery(n,o,p,q,r,s,t,u,v,w,x) = new(n,o,p,q,r,s,t,u,v,w,x,0.0,24)
    battery(n,o,p,q,r,s,t,u,v,w) = new(n,o,p,q,r,s,t,u,v,w,1.0,0.0,24)
    battery(n,o,p,q,r,s,t,u) = new(n,o,p,q,r,s,t,u,1.0,1.0,1.0,0.0,24)
    battery(n,o,p,q,r,s,t) = new(n,o,p,q,r,s,t,"new",1.0,1.0,1.0,0.0,24)
   
    # Checks
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = (!(0.0 <= y <= 1.0) || (z < 0)) ? error("FOR and/or MTTR values passed are not allowed") : new(n,o,p,q,r,s,t,u,v,w,x,y,z)
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = !(u in ["old","new"]) ? error("Unidentified legacy passed") : new(n,o,p,q,r,s,t,u,v,w,x,y,z)
    battery(n,o,p,q,r,s,t,u,v,w,x,y,z) = !(0 <= o<= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(n,o,p,q,r,s,t,u,v,w,x,y,z)
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
    gen_storage(name = "stor_1", N = 10, rg = "reg_1",type = "4-hour", c_cap = fill(10.0,1,N), dis_cap = fill(10.0,1,N), energy_cap = fill(40.0,1,N), infl = fill(10.0,1,N), 
                g_with_cap = fill(10.0,1,N),g_inj_cap = fill(10.0,1,N), lg = "new", chr_eff = 0.9, dis_eff = 1.0, cry_eff = 1.0, f_or = 0.0, mttr = 24) = 
                new(name,N,rg,type, c_cap,dis_cap,energy_cap, infl, g_with_cap, g_inj_cap, lg, chr_eff, dis_eff,cry_eff, f_or, mttr)

    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x) = new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,0.0,24)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w) = new(k,l,m,n,o,p,q,r,s,t,u,v,w,1.0,0.0,24)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u) = new(k,l,m,n,o,p,q,r,s,t,u,1.0,1.0,1.0,0.0,24)
    gen_storage(k,l,m,n,o,p,q,r,s,t) = new(k,l,m,n,o,p,q,r,s,t,"new",1.0,1.0,1.0,0.0,24)
   
    # Checks
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = (!(0.0 <= y <= 1.0) || (z<0)) ? error("FOR and/or MTTR values passed are not allowed") : new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(u in ["old","new"]) ? error("Unidentified legacy passed") : new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
    gen_storage(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z) = !(0 <= l <= 8784) ? error("Check the PRAS timesteps (N) passed.") : new(k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z)
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
    return getfield(outage_to_rate((stor.FOR,stor.MTTR)),:λ)
end

function get_μ(stor::STOR) where {STOR <: storage}
    return getfield(outage_to_rate((stor.FOR,stor.MTTR)),:μ)
end

function get_category(stor::STOR) where {STOR <: storage}
    return stor.legacy*stor.type
end

function get_charge_capacity(stor::battery)
    return fill(round(Int,stor.charge_cap),1,stor.N)
end

function get_charge_capacity(stor::gen_storage)
    return round.(Int,stor.charge_cap)
end

function get_discharge_capacity(stor::battery)
    return fill(round(Int,stor.discharge_cap),1,stor.N)
end

function get_discharge_capacity(stor::gen_storage)
    return round.(Int,stor.discharge_cap)
end

function get_energy_capacity(stor::battery)
    return fill(round(Int,stor.energy_cap),1,stor.N)
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
    leg_stor_idxs = findall(getfield.(stors,:legacy) .== leg)
    if isnothing(eg_stor_idxs)
        @warn "No generators with this legacy"
    else
        return stors[leg_stor_idxs]
    end
end
