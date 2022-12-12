function clean_names(input_vec::Vector)
    for (idx,a) in enumerate(input_vec)
        if occursin("*",a)
            input_vec[idx] = match(r"\*([a-zA-Z]+-*[a-zA-Z]*)_*", a)[1]
        end
    end 
    return input_vec
end

function expand_types(input_vec::Vector,N::Int64)
    add_names = [];
    for a in input_vec
        for n in 1:N
            mystr = string(a*"_"*string(n))
            push!(add_names,mystr)
        end
    end
    return vcat(input_vec,add_names)
end

function sort_gens(gen_regionnames::Vector,regionnames::Vector,n_regions::Int)
    regionlookup = Dict(n=>i for (i, n) in enumerate(regionnames))
    gen_regions = getindex.(Ref(regionlookup), gen_regionnames)
    region_order = sortperm(gen_regions)
    region_gen_idxs = makeidxlist(gen_regions[region_order], n_regions)
    return (region_gen_idxs,region_order)
end

#borrowed from Gord's PRAS Utils fxn
#I guess I could just call this straight from PRAS, but let's leave it for now
function makeidxlist(collectionidxs::Vector{Int}, n_collections::Int)

    n_assets = length(collectionidxs)

    idxlist = Vector{UnitRange{Int}}(undef, n_collections)
    active_collection = 1
    start_idx = 1
    a = 1

    while a <= n_assets
       if collectionidxs[a] > active_collection
            idxlist[active_collection] = start_idx:(a-1)       
            active_collection += 1
            start_idx = a
       else
           a += 1
       end
    end

    idxlist[active_collection] = start_idx:n_assets       
    active_collection += 1

    while active_collection <= n_collections
        idxlist[active_collection] = (n_assets+1):n_assets
        active_collection += 1
    end

    return idxlist

end

#Borrowed & modified from PLEXOS2PRAS 
#only have FOR right now, so just make simple assumption about MTTR being 24 hours (IDK, whatever)
function FOR_to_transitionprobs(for_raw::Float64) 

    # raw MTTR is in hours, raw FOR is a fraction
    mttr = 24. #input for now 
    μ = 1 ./ mttr
    # μ[mttrs .== 0] .= one(V) # Interpret zero MTTR as μ = 1.
    # fors = for_raw ./ 100
    λ = μ .* for_raw ./ (1 .- for_raw)
    # λ[fors .== 0] .= zero(V) # Interpret zero FOR as λ = 0.
    return λ, μ

end

function run_pras_system(sys::PRAS.SystemModel,sample::Int)
    shortfalls,flows = PRAS.assess(sys,PRAS.SequentialMonteCarlo(samples=sample),PRAS.Shortfall(),PRAS.Flow())
    println(PRAS.LOLE(shortfalls))
    return shortfalls,flows
end

### disaggregation of capacity is not yet implemented ###

function disagg_existing_capacity(eia_df::DataFrames.DataFrame,built_capacity::Int,tech::String,pca::String,Year::Int)
    tech_ba_year_existing = eia_df[(eia_df.tech.==tech) .& (eia_df.reeds_ba.==pca) .& (eia_df.RetireYear.>=Year) .& (eia_df.StartYear.<=Year), :];
    remaining_capacity = built_capacity;
    existing_gens = [];
    tech_len = length(tech_ba_year_existing[!,"cap"]);
    max_cap = maximum(eia_df[!,"cap"])
    avg_cap = mean(eia_df[!,"cap"])

    for (idx,built_cap) in enumerate(tech_ba_year_existing[!,"cap"])
        int_built_cap = floor.(Int,built_cap);
        if int_built_cap < remaining_capacity
            remaining_capacity = remaining_capacity - int_built_cap;
            gen_name = tech*"_"*pca*"_"*string.(idx);
            #create a generator using types.jl
            existing_gen = thermal_gen(gen_name,8760,pca,int_built_cap,tech,"Existing",.05,24) #here we'd actually want to make a choice about which struct to send the generator to, and pass N and FOR
            push!(existing_gens,existing_gen);
        else
            @info "on $tech $idx of $tech_len in $pca, $remaining_capacity MW remains but can't be built. Build as one last existing generator of $remaining_capacity MW for now"
            gen_name = tech*"_"*pca*"_"*string.(idx);
            existing_gen = thermal_gen(gen_name,8760,pca,remaining_capacity,tech,"Existing",.05,24);
            remaining_capacity = 0;
            break
        end
        
    end
    #whatever remains, we want to build as new capacity
    @info "overall, for $tech $pca, $remaining_capacity of $built_capacity MW will be new build with target average size $avg_cap MW"
    if remaining_capacity > 0
        new_gens = disagg_new_capacity(remaining_capacity,floor.(Int,avg_cap),floor.(Int,max_cap),tech,pca,Year);
    else
        new_gens = [];    
    
    end
    #pass each region-gentype capacity, filter against EIA-NEMS and retirement info
    #then return disaggregated region-gentype existing capacity and newbuild remainder
    #push remaining capacity to new capacity fxn
    
    return (existing_gens,new_gens)
    # return tech_ba_year_existing
end

function disagg_new_capacity(new_capacity::Int,avg::Int,max::Int,tech::String,pca::String,Year::Int)
    # cap_out = [];
    gens = floor.(Int,new_capacity/avg);
    remainder = new_capacity-(gens*avg);
    addtl_cap_per_gen = floor.(Int,remainder/gens);
    per_gen_cap = avg+addtl_cap_per_gen;
    small_remainder = new_capacity-(gens*per_gen_cap)
    @info "after remainder peanut butter, average capacity is now $per_gen_cap MW for $gens new generators, with an additional $small_remainder MW unit to be built"
    cap_out = [thermal_gen(tech*"_"*pca*"_new_"*string.(i),8760,pca,per_gen_cap,tech,"New",.05,24) for i in range(1,gens)]; #then make the gens
    push!(cap_out,thermal_gen(tech*"_"*pca*"_new_"*string.(gens+1),8760,pca,small_remainder,tech,"New",.05,24)); #integer remainder is made into a tiny gen
    return cap_out
end