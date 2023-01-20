#####################################################
# Luke
# NREL
# November 2022
# ReEDS2PRAS for NTP
# Make a PRAS System from ReEDS H5s and CSVs

function reeds_to_pras(ReEDSfilepath::String, year::Int64, NEMS_path::String, N::Int, WEATHERYEAR::Int)
    #######################################################
    # Loading the necessary mapping files and data
    #######################################################
    min_year = 2007 #for now, read in later
    if WEATHERYEAR ∉ [min_year,2008,2009,2010,2011,2012,2013] 
        error("The weather year $WEATHERYEAR is not a valid VG profile year. Should be an Int in $min_year-$(min_year+6) currently")
    end
    ReEDS_data = ReEDSdata(ReEDSfilepath,year)
    
    #######################################################
    # Create objects for a PRAS system
    #######################################################
    all_lines,regions,all_gens,storages_array,genstor_array = create_objects(ReEDS_data,WEATHERYEAR,N,year,NEMS_path,min_year)

    @info "...objects are created, writing to PRAS system"
    return create_pras_system(regions,all_lines,all_gens,storages_array,genstor_array,N,year)
end

function create_objects(ReEDS_data,WEATHERYEAR::Int,N::Int,year::Int,NEMS_path::String,min_year::Int)
    
    #######################################################
    # Create regions and associate load profiles
    #######################################################
    @info "Processing regions..."
    region_array = process_regions_and_load(ReEDS_data,WEATHERYEAR,N)
    
    #######################################################
    # Create line objects & add VSC-related regions, if applicable
    #######################################################
    @info "Processing lines..."
    lines = process_lines(ReEDS_data,get_name.(region_array),year,N)
    lines,regions = process_vsc_lines(lines,region_array)
    
    #######################################################
    # Create Generator Objects
    # **TODO: Should 0 MW generators be allowed after disaggregation?
    # **TODO: Should hydro be split out as a generator-storage?
    # **TODO: is it important to also handle planned outages?
    #######################################################
    @info "splitting thermal, storage, vg generator types from installed ReEDS capacities..."
    thermal,storage,variable_gens = split_generator_types(ReEDS_data,year)

    @info "reading in ReEDS generator-type forced outage data..."
    forced_outage_data = get_forced_outage_data(ReEDS_data)
    forced_outage_dict = Dict(forced_outage_data[!,"ResourceType"] .=> forced_outage_data[!,"FOR"])

    @info "reading thermals..."
    thermal_gens = process_thermals_with_disaggregation(thermal,forced_outage_dict,N,year,NEMS_path)
    @info "reading vg..."
    gens_array = process_vg(thermal_gens,variable_gens,forced_outage_dict,ReEDS_data,year,WEATHERYEAR,N,min_year)

    #######################################################
    # Create Storage Objects
    #######################################################
    @info "Processing Storages..."
    storage_array = process_storages(storage,forced_outage_dict,ReEDS_data,N,get_name.(regions),year)

    #######################################################
    # Creating GenStor Objects
    #######################################################
    @info "Processing GeneratorStorages [EMPTY FOR NOW].."
    genstor_array = process_genstors(get_name.(regions),N)

    return lines,regions,gens_array,storage_array,genstor_array
end

function create_pras_system(regions::Vector{Region},lines::Vector{Line},gens::Vector{<:Generator},storages::Vector{<:Storage},gen_stors::Vector{<:Gen_Storage},N::Int,year::Int)

    first_ts = TimeZones.ZonedDateTime(year, 01, 01, 00, TimeZones.tz"EST") #switched to US/EST from UTC
    last_ts = first_ts + Dates.Hour(N-1)
    my_timestamps = StepRange(first_ts, Dates.Hour(1), last_ts)

    sorted_lines,interface_reg_idxs,interface_line_idxs = get_sorted_lines(lines,regions)
    pras_lines,pras_interfaces = make_pras_interfaces(sorted_lines,interface_reg_idxs,interface_line_idxs,regions)
    pras_regions = PRAS.Regions{N,PRAS.MW}(get_name.(regions), reduce(vcat,get_load.(regions)))
    ##
    sorted_gens, gen_idxs = get_sorted_components(gens,get_name.(regions))
    capacity_matrix = reduce(vcat,get_capacity.(sorted_gens))
    λ_matrix = reduce(vcat,get_λ.(sorted_gens))
    μ_matrix = reduce(vcat,get_μ.(sorted_gens))
    pras_gens = PRAS.Generators{N,1,PRAS.Hour,PRAS.MW}(get_name.(sorted_gens),get_type.(sorted_gens),capacity_matrix,λ_matrix,μ_matrix)
    ##
    storages, stor_idxs = get_sorted_components(storages,regions);

    stor_charge_cap_array = reduce(vcat,get_charge_capacity.(storages))
    stor_discharge_cap_array = reduce(vcat,get_discharge_capacity.(storages))
    stor_energy_cap_array = reduce(vcat,get_energy_capacity.(storages))
    stor_chrg_eff_array = reduce(vcat,get_charge_efficiency.(storages))
    stor_dischrg_eff_array = reduce(vcat,get_discharge_efficiency.(storages))
    stor_cryovr_eff = reduce(vcat,get_carryover_efficiency.(storages))
    λ_stor = reduce(vcat,get_λ.(storages))
    μ_stor = reduce(vcat,get_μ.(storages))
    pras_storages = PRAS.Storages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(get_name.(storages),get_type.(storages),
                                                stor_charge_cap_array,stor_discharge_cap_array,stor_energy_cap_array,
                                                stor_chrg_eff_array,stor_dischrg_eff_array, stor_cryovr_eff,
                                                λ_stor,μ_stor)
    ##
    sorted_gen_stors, genstor_idxs  = get_sorted_components(gen_stors,regions)

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

    gen_stors = PRAS.GeneratorStorages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(gen_stor_names,gen_stor_cats,gen_stor_cap_array, gen_stor_dis_cap_array, gen_stor_enrgy_cap_array,
                                                                           gen_stor_chrg_eff_array, gen_stor_dischrg_eff_array, gen_stor_carryovr_eff_array,gen_stor_inflow_array,
                                                                           gen_stor_grid_withdrawl_array, gen_stor_grid_inj_array,gen_stor_λ,gen_stor_μ)

    return PRAS.SystemModel(pras_regions, pras_interfaces, pras_gens, gen_idxs, pras_storages, stor_idxs, gen_stors, genstor_idxs,
                            pras_lines,interface_line_idxs,my_timestamps)
                                        
end

function process_regions_and_load(ReEDS_data,weather_year::Int, N::Int)

    load_info = get_load_file(ReEDS_data)
    load_data = load_info["block0_values"]
    regions = load_info["block0_items"]
    #To-Do: can get a bug/error here if ReEDS case lacks multiple load years
    slicer = findfirst(isequal(weather_year), load_info["axis1_level1"]); #slices based on input weather year
    indices = findall(i->(i==(slicer-1)),load_info["axis1_label1"]); #I think b/c julia indexes from 1 we need -1 here
    
    load_year = load_data[:,indices[1:N]]; #should be regionsX8760, which is now just enforced

    return [Region(r,N,floor.(Int,load_year[idx,:])) for (idx,r) in enumerate(regions)]
end

function process_lines(ReEDS_data,regions::Vector{<:AbstractString},year::Int,N::Int)

    line_base_cap_data = get_line_capacity_data(ReEDS_data)
    converter_capacity_data = get_converter_capacity_data(ReEDS_data)
    converter_capacity_dict = Dict(converter_capacity_data[!,"r"] .=> converter_capacity_data[!,"MW"])

    #add 0 converter capacity for regions that lack a converter
    for reg in regions
        if !(reg in keys(converter_capacity_dict))
            @info "$reg does not have VSC converter capacity, so adding a 0"
            converter_capacity_dict[string(reg)] = 0
        end
    end

    system_line_idx = []
    region_idxs = Dict(regions .=> range(1,length(regions)))
    for (idx,row) in enumerate(eachrow(line_base_cap_data))
        from_idx = region_idxs[row.r]
        to_idx = region_idxs[row.rr]
        if (~(isnothing(from_idx)) && ~(isnothing(to_idx)) && (from_idx < to_idx))
            push!(system_line_idx,idx)
        end
    end
    #order is assumed preserved in splitting these dfs for now but should likely be checked
    system_line_naming_data = line_base_cap_data[system_line_idx,:] # all line capacities

    ###TODO: THERE IS A BETTER WAY HERE BUT I NEED TO UNDERSTAND SUBSET IN JULIA 
    # function keep_line(from_pca, to_pca, region_idxs)
    #     from_idx = region_idxs[from_pca]
    #     to_idx = region_idxs[to_pca]
    #     return ~isnothing(from_idx) && ~isnothing(to_idx) && (from_idx < to_idx)
    # end
    # system_line_naming_data = DataFrames.subset(line_base_cap_data,keep_line(line_base_cap_data.r, line_base_cap_data.rr, region_idxs))

    # closure - keep_line() function has region line information embedded and instructions. 

    system_line_naming_data = DataFrames.combine(DataFrames.groupby(system_line_naming_data, ["r","rr","trtype"]), :MW => sum) #split-apply-combine b/c some lines have same name convention

    lines_array = Line[]
    for row in eachrow(system_line_naming_data)
        forward_cap = sum(line_base_cap_data[(line_base_cap_data.r.==row.r) .& (line_base_cap_data.rr.==row.rr) .& (line_base_cap_data.trtype.==row.trtype),"MW"])
        backward_cap = sum(line_base_cap_data[(line_base_cap_data.r.==row.rr) .& (line_base_cap_data.rr.==row.r) .& (line_base_cap_data.trtype.==row.trtype),"MW"])

        name = "$(row.r)_$(row.rr)_$(row.trtype)"        
        @info "a line $name, with $forward_cap MW fwd and $backward_cap bckwrd in $(row.trtype)"
        if row.trtype!="VSC"
            push!(lines_array,Line(name,N,row.trtype,row.r,row.rr,forward_cap,backward_cap,"Existing",0.0,24))#for and mttr will be defaults
        else
            push!(lines_array,Line(name,N,row.trtype,row.r,row.rr,forward_cap,backward_cap,"Existing",0.0,24,true,Dict(row.r => converter_capacity_dict[string(row.r)], row.rr => converter_capacity_dict[string(row.rr)])))
        end
    end
    return lines_array
end

function expand_vg_types!(vgl::Vector{<:AbstractString},vgt::Vector)
    #TODO: vgt/vgl check
    for l in vgl
        @assert occursin(l,join(vgt))
    end 
    return vec(["$(a)" for a in vgt])
end

function split_generator_types(ReEDS_data,year::Int64)

    tech_types_data = get_technology_types(ReEDS_data)
    capacity_data = get_ICAP_data(ReEDS_data)
    vg_resource_types = get_valid_resources(ReEDS_data)

    vg_types = DataFrames.dropmissing(tech_types_data,:VRE)[:,"Column1"]
    vg_types = vg_types[vg_types .!= "csp-ns"]#csp-ns causes problems, so delete for now
    storage_types = DataFrames.dropmissing(tech_types_data,:STORAGE)[:,"Column1"]

    #clean vg/storage capacity on a regex, though there might be a better way...    
    clean_names!(vg_types)
    clean_names!(storage_types)
    #vg_types = expand_types!(vg_types,15) #expand so names will match
    vg_types = expand_vg_types!(vg_types,unique(vg_resource_types.i)) 

    vg_capacity = capacity_data[(findall(in(vg_types),capacity_data.i)),:]
    storage_capacity = capacity_data[(findall(in(storage_types),capacity_data.i)),:]
    thermal_capacity = capacity_data[(findall(!in(vcat(vg_types,storage_types)),capacity_data.i)),:]
    return(thermal_capacity,storage_capacity,vg_capacity)
end

function process_thermals_with_disaggregation(thermal_builds::DataFrames.DataFrame,FOR_dict::Dict,N::Int,year::Int,NEMS_path::String)#FOR_data::DataFrames.DataFrame,
    thermal_builds = thermal_builds[(thermal_builds.i.!= "csp-ns"), :] #csp-ns is not a thermal; just drop in for rn
    thermal_builds = DataFrames.combine(DataFrames.groupby(thermal_builds, ["i","r"]), :MW => sum) #split-apply-combine to handle differently vintaged entries
    EIA_db = Load_EIA_NEMS_DB(NEMS_path)
    
    all_generators = Generator[]
    lowercase_FOR_dict = Dict([lowercase(i) for i in keys(FOR_dict)] .=> values(FOR_dict))
    for row in eachrow(thermal_builds) #this loop gets the FOR for each build/tech
        # @info "$(row.i) $(row.r) $(row.MW_sum) MW to disaggregate..."
        if row.i in keys(FOR_dict)#native_FOR_data
            gen_for = FOR_dict[row.i]
        elseif row.i in keys(lowercase_FOR_dict)
            gen_for = lowercase_FOR_dict[row.i]
        else
            gen_for = 0.05
            @info "for $(row.i) $(row.r), no gen_for is found in data, so $gen_for is used"
        end
        generator_array = disagg_existing_capacity(EIA_db,floor(Int,row.MW_sum),string(row.i),string(row.r),gen_for,N,year)
        append!(all_generators,generator_array)
    end
    return all_generators
end

function process_vg(generators_array::Vector{<:ReEDS2PRAS.Generator},vg_builds::DataFrames.DataFrame,FOR_dict::Dict,ReEDS_data,year::Int,weather_year::Int,N::Int,min_year::Int)
    #data loads
    region_mapper_df = get_region_mapping(ReEDS_data)
    region_mapper_dict = Dict(region_mapper_df[!,"rs"] .=> region_mapper_df[!,"*r"])
    cf_info = get_vg_cf_data(ReEDS_data)#load is now picked up from augur
    
    vg_profiles = cf_info["block0_values"]
    start_idx = (weather_year-min_year)*N

    vg_builds = DataFrames.combine(DataFrames.groupby(vg_builds, ["i","r"]), :MW => sum) #split-apply-combine to handle differently vintaged entries

    for row in eachrow(vg_builds)
        category = string(row.i)
        name = "$(category)_$(string(row.r))"
        
        if string(row.r) in keys(region_mapper_dict)
            region = region_mapper_dict[string(row.r)]
        else
            region = string(row.r)
        end

        profile_index = findfirst.(isequal.(name), (cf_info["axis0"],))[1]
        size_comp = size(vg_profiles)[1]
        profile = vg_profiles[profile_index,(start_idx+1):(start_idx+N)]
        if category in keys(FOR_dict)
            gen_for = FOR_dict[category]
        else
            gen_for = .05
        end
        name = "$(name)_"
        
        push!(generators_array,VG_Gen(name,N,region,maximum(profile),profile,category,"New",gen_for,24)) 
    end
    return generators_array
end

function process_storages(storage_builds::DataFrames.DataFrame,FOR_dict::Dict,ReEDS_data,N::Int,regions::Vector{<:AbstractString},year::Int64)
    
    storage_energy_capacity_data = get_storage_energy_capacity_data(ReEDS_data)
    energy_capacity_df = DataFrames.combine(DataFrames.groupby(storage_energy_capacity_data, ["i","r"]), :MWh => sum) #split-apply-combine to handle differently vintaged entries

    storages_array = Storage[];
    for (idx,row) in enumerate(eachrow(storage_builds))
        name = "$(string(row.i))_$(string(row.r))"
        if string(row.i) in keys(FOR_dict)
            gen_for = FOR_dict[string(row.i)]
        else
            gen_for = 0.0;
            @info "did not find FOR for storage $name $(row.r) $(row.i), so setting FOR to default value $gen_for"
        end 
        name = "$(name)_"#append for later matching
        push!(storages_array,Battery(name,N,string(row.r),string(row.i),row.MW,row.MW,energy_capacity_df[idx,"MWh_sum"],"New",1,1,1,gen_for,24)) 
    end
    return storages_array
end

function process_genstors(regions::Vector{<:AbstractString},N::Int)
    
    gen_stors = Gen_Storage[Gen_Storage("blank_1", N, regions[1], "blank_genstor", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)] #empty for now
    
    return gen_stors
end