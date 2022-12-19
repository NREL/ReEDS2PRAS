function capacity_checker(capacity_data::DataFrames.DataFrame,gentype::String,region::String)
    #@assert gentype in unique(capacity_data.i) #check valid pass
    #@assert region in unique(capacity_data.r) #check valid pass
    capacity_data_subset = capacity_data[(capacity_data.i.==gentype) .& (capacity_data.r.==region), :];
    return sum(capacity_data_subset[!,"MW"])# sum the capacity for the input region/gentype
end

function PRAS_generator_capacity_checker(pras_system,gentype::String,region::String)
    name_vec = occursin.(gentype,pras_system.generators.names);
    reg_vec = occursin.(region,pras_system.generators.names);
    out_vec = .*(name_vec,reg_vec); 
    retained_gens = [];
    for (idx,val) in enumerate(out_vec) #there has to be a better, way, but for now
        if val==1
            push!(retained_gens,idx)
        end
    end
    return sum(psys.generators.capacity[Int.(retained_gens)])
end

function PRAS_storage_capacity_checker(pras_system,gentype::String,region::String)
    name_vec = occursin.(gentype,pras_system.storages.names);
    reg_vec = occursin.(region,pras_system.storages.names);
    out_vec = .*(name_vec,reg_vec); 
    retained_gens = [];
    for (idx,val) in enumerate(out_vec) #there has to be a better, way, but for now
        if val==1
            push!(retained_gens,idx)
        end
    end
    return sum(psys.storages.discharge_capacity[Int.(retained_gens)])
end

function compare_generator_capacities(psys,ReEDSfilepath,Year)
    #first actually have to load in case-level capacity data, which may not be passed
    ReEDS_data = ReEDS2PRAS.ReEDSdata(ReEDSfilepath,Year);
    tech_types_data = ReEDS2PRAS.get_technology_types(ReEDS_data);
    techs = tech_types_data[!,"Column1"];
    capacity_data = ReEDS2PRAS.get_ICAP_data(ReEDS_data);
    for gentype in techs
        for region in psys.regions.names
            # @info "comparing $gentype-$region input/output capacities..."
            v1 = capacity_checker(capacity_data,gentype,region);
            v2a = PRAS_generator_capacity_checker(psys,gentype,region);
            v2b = PRAS_storage_capacity_checker(psys,gentype,region);
            v2 = v2a+v2b;
            if v1 != 0 || v2 != 0
                @info "for $gentype-$region input is $v1 MW, output is $v2 MW"
            end
        end
    end
end