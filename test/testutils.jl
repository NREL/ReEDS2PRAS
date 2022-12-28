function capacity_checker(capacity_data::DataFrames.DataFrame,region_map::DataFrames.DataFrame,gentype::String,region::String)
    #@assert gentype in unique(capacity_data.i) #check valid pass
    #@assert region in unique(capacity_data.r) #check valid pass
    for idx in range(1,DataFrames.nrow(capacity_data))
        slicer = findfirst(isequal(string(capacity_data[idx,"r"])), region_map[:,"rs"]);
        if !isnothing(slicer)
            capacity_data[idx,"r"] = region_map[slicer,"*r"];
        end
    end
    
    capacity_data_subset = capacity_data[(capacity_data.i.==gentype) .& (capacity_data.r.==region), :];
    return sum(capacity_data_subset[!,"MW"])# sum the capacity for the input region/gentype
end

function PRAS_generator_capacity_checker(pras_system,gentype::String,region::String)

    name_vec = -abs.(cmp.(gentype,pras_system.generators.categories)).+1; #exact match is needed to exclude ccs
    
    reg_vec = occursin.(region*"_",pras_system.generators.names);
    out_vec = .*(name_vec,reg_vec); 

    retained_gens = [];
    for (idx,val) in enumerate(out_vec) #there has to be a better, way, but for now
        if val==1
            push!(retained_gens,idx)
        end
    end

    return sum(maximum.(eachrow(pras_system.generators.capacity[Int.(retained_gens),:])))
end

function PRAS_storage_capacity_checker(pras_system,gentype::String,region::String)
    name_vec = occursin.(gentype,pras_system.storages.names);
    reg_vec = occursin.(region*"_",pras_system.storages.names);
    out_vec = .*(name_vec,reg_vec); 
    retained_gens = [];
    for (idx,val) in enumerate(out_vec) #there has to be a better, way, but for now
        if val==1
            push!(retained_gens,idx)
        end
    end
    return sum(maximum.(eachrow(pras_system.storages.discharge_capacity[Int.(retained_gens),:])))
end

function clean_gentype(input_name::String)
    if occursin("*",input_name)
        input_name = string(match(r"\*([a-zA-Z]+-*[a-zA-Z]*)_*", input_name)[1]);
        input_vec = expand_vg_types([input_name],15);
    else
        input_vec = [input_name];
    end
    return input_vec
end

function expand_vg_types(input_vec::Vector,N::Int64)
    add_names = [];
    for a in input_vec
        for n in 1:N
            mystr = string(a*"_"*string(n));
            push!(add_names,mystr);
        end
    end
    return vcat(input_vec,add_names)
end

function compare_generator_capacities(pras_system,ReEDSfilepath,Year)
    #first actually have to load in case-level capacity data, which may not be passed
    ReEDS_data = ReEDS2PRAS.ReEDSdata(ReEDSfilepath,Year);

    capacity_data = ReEDS2PRAS.get_ICAP_data(ReEDS_data);
    region_mapper_df = ReEDS2PRAS.get_region_mapping(ReEDS_data);

    for gentype in unique(capacity_data.i)
        gentype = string(gentype);
        for region in pras_system.regions.names
            v1 = capacity_checker(capacity_data,region_mapper_df,gentype,region); #need to split out numbering for vg...
            v2a = PRAS_generator_capacity_checker(pras_system,gentype,region);
            v2b = PRAS_storage_capacity_checker(pras_system,gentype,region);
            v2 = v2a+v2b;
            if v1 != 0 || v2 != 0
                if abs(v1-v2)>1
                    @info "for $gentype-$region input is $v1 MW, output is $v2 MW. Mismatch!!"
                end
            end
        end
    end
end

function compare_line_capacities(pras_system,ReEDSfilepath,Year)
    ReEDS_data = ReEDS2PRAS.ReEDSdata(ReEDSfilepath,Year);
    line_df = ReEDS2PRAS.get_line_capacity_data(ReEDS_data);

    for (r1,r2,type,MW) in zip(line_df[!,"r"],line_df[!,"rr"],line_df[!,"trtype"],line_df[!,"MW"])
        r1_vec = occursin.(r1*"_",pras_system.lines.names);
        r2_vec = occursin.(r2*"_",pras_system.lines.names);
        type_vec = occursin.(type,pras_system.lines.names);
        out_vec = .*(r1_vec,r2_vec,type_vec); 
        retained_lines = [];
        mw_sum = 0;
        for (idx,val) in enumerate(out_vec) #there has to be a better, way, but for now
            if val==1
                push!(retained_lines,idx)
                mw_sum = mw_sum+MW;
            end
        end
        mw_out_fwd,mw_out_bck = get_pras_line_capacity(pras_system,Int.(retained_lines))
        # if parse(Int64, r1[2:end])<parse(Int64, r2[2:end]) #this is a fwd line cap
        if r1 < r2 #actually want to compare strings to get proper ordering
            @info "forward in MW is $mw_sum out forward is $mw_out_fwd, for $r1 to $r2, $type"
            @assert abs(mw_sum-mw_out_fwd)<=1
        else #bckwd line cap
            #mw_out_fwd,mw_out_bck = get_pras_line_capacity(pras_system,Int.(retained_lines))
            @info "backward in MW is $mw_sum out backward is $mw_out_bck for $r2 to $r1, $type"
            @assert abs(mw_sum-mw_out_bck)<=1
        end
    end
end

function get_pras_line_capacity(pras_system,idx_list::Vector)
    #TODO: check backward capacity? 
    fwd = sum(maximum.(eachrow(pras_system.lines.forward_capacity[idx_list,:]))); #get pras line fwd capacity
    bck = sum(maximum.(eachrow(pras_system.lines.backward_capacity[idx_list,:]))); #get pras line fwd capacity
    return fwd,bck
end