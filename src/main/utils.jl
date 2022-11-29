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