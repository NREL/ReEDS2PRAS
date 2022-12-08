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

    if (for_gen >1.0)
        for_gen = for_gen/100
    end

    if (mttr != 0)
        μ = 1 / mttr
    else
        μ = 0.0
    end
    λ = (μ * for_gen) / (1 - for_gen)
    #λ = for_gen

    return (λ = λ, μ = μ)
end
# Regions
struct region
    name::String
    N::Int64
    load::Vector{Float64}

    region(name = "region_1", N = 10, load = ones(Float64,N)) = new(name,N,load)
    region(name = "region_1", N =10) = new(name,N, zeros(Float64,N))

    region(x,y,z,) = (length(z) !== y) ? error("The length of the region load should be equal to PRAS timesteps (N)") : new(x,y,z)
end

function get_name(reg::region)
    return reg.name
end

function get_load(reg::region)
    return reg.load
end

# Generators
struct generator
    name::String
    FOR::Float64
    MTTR:: Float64

    generator(name = "gen_1", f_or = 0.1, mttr = 24) = new(name,f_or,mttr)
    generator(name = "gen_1", f_or = 0.1) = new(name,f_or,24)
    generator(name = "gen_1") = new(name,0.0,24)

    generator(x,y,z,) = (!(0.0<=y<=1.0) || (z<0)) ? error("FOR and/or MTTR values passed are not allowed") : new(x,y,z)
end

function get_name(gen::GEN) where {GEN <: generator}
    gen.name
end