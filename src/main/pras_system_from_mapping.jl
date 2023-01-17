#####################################################
# Luke
# NREL
# November 2022
# ReEDS2PRAS for NTP
# Make a PRAS System from ReEDS H5s and CSVs

function reeds_to_pras(ReEDSfilepath::String, year::Int64, NEMS_path::String)
    #######################################################
    # Loading the necessary mapping files and data
    #######################################################
    #check validity of input weather and ReEDS year
    WEATHERYEAR = 2012 #just for now, force this
    if WEATHERYEAR ∉ [2007,2008,2009,2010,2011,2012,2013] 
        error("The weather year $WEATHERYEAR is not a valid VG profile year. Should be an Int in 2007-2013 currently")
    end
    ReEDS_data = ReEDSdata(ReEDSfilepath,year);

    #######################################################
    # Load the mapping metadata JSON file
    # TODO: depend on PRAS regions, not region_array in later functions
    #######################################################
    
    VSC_append_str,N,region_array,new_regions = regions_and_load(ReEDS_data,WEATHERYEAR)

    #######################################################
    # PRAS lines
    #######################################################

    new_lines,new_interfaces,interface_line_idxs = process_lines(ReEDS_data,get_name.(region_array),year,N,VSC_append_str);
    
    #######################################################
    # PRAS Generators
    # **TODO: Should 0 MW generators be allowed after disaggregation?
    # **TODO: Should hydro be split out as a generator-storage?
    # **TODO: is it important to also handle planned outages?
    #######################################################

    @info "splitting thermal, storage, vg generator types from installed ReEDS capacities..."
    thermal,storage,vg = split_generator_types(ReEDS_data,year);

    @info "reading in ReEDS generator-type forced outage data..."
    forced_outage_data = get_forced_outage_data(ReEDS_data);

    @info "reading thermals..."
    therm_gens = process_thermals_with_disaggregation(thermal,forced_outage_data,N,year,NEMS_path);
    @info "reading vg..."
    all_gens = process_vg(therm_gens,vg,forced_outage_data,ReEDS_data,year,WEATHERYEAR,N);
    @info "...vg read, translating to PRAS gens"
    new_generators,area_gen_idxs = all_gens_to_pras(all_gens,region_array,N)

    #######################################################
    # PRAS Timestamps
    #######################################################
    @info "Processing PRAS timestamps..."
    first_ts = TimeZones.ZonedDateTime(year, 01, 01, 00, TimeZones.tz"UTC"); #US/Eastern switch to EST/EDT
    last_ts = first_ts + Dates.Hour(N-1);
    my_timestamps = StepRange(first_ts, Dates.Hour(1), last_ts);

    #######################################################
    # PRAS Storages
    #######################################################
    @info "Processing Storages..."
    area_stor_idxs,new_storage = process_storages(storage,forced_outage_data,ReEDS_data,N,get_name.(region_array),year)
    
    #######################################################
    # PRAS GeneratorStorages
    #######################################################
    @info "Processing GeneratorStorages [EMPTY FOR NOW].."
    new_gen_stors,area_genstor_idxs = process_genstors(get_name.(region_array),N)

    pras_system = PRAS.SystemModel(new_regions, new_interfaces, new_generators, area_gen_idxs, new_storage, area_stor_idxs, new_gen_stors,
                                    area_genstor_idxs, new_lines,interface_line_idxs,my_timestamps);
    return pras_system
end

function all_gens_to_pras(all_gens::Vector{<:ReEDS2PRAS.Generator},region_array::Vector,N::Int)
    sorted_gens, area_gen_idxs = get_sorted_components(all_gens,get_name.(region_array)); #TODO: is typing still wrong?
    
    capacity_matrix = reduce(vcat,get_capacity.(sorted_gens));
    λ_matrix = reduce(vcat,get_λ.(sorted_gens));
    μ_matrix = reduce(vcat,get_μ.(sorted_gens));

    return (PRAS.Generators{N,1,PRAS.Hour,PRAS.MW}(get_name.(sorted_gens),get_type.(sorted_gens),capacity_matrix,λ_matrix,μ_matrix),area_gen_idxs)
end

function regions_and_load(ReEDS_data::CEMdata,WEATHERYEAR::Int)
    @info "Fetching ReEDS case data to build PRAS System..."

    load_info = get_load_file(ReEDS_data);
    load_data = load_info["block0_values"];
    regions = load_info["block0_items"];
    #To-Do: can get a bug/error here if ReEDS case lacks multiple load years
    slicer = findfirst(isequal(WEATHERYEAR), load_info["axis1_level1"]); #2012 weather year, Gbut, in general, this should be a user input
    indices = findall(i->(i==(slicer-1)),load_info["axis1_label1"]); #I think b/c julia indexes from 1 we need -1 here
    
    load_year = load_data[:,indices[1:8760]]; #should be regionsX8760, which is now just enforced
    
    #should the regions be processed from the generators?
    @info "Processing Areas in the mapping file into PRAS regions..."
    num_areas = length(regions); 
    N = size(load_year)[2];#DataFrames.DataFrames.nrow(load_data)
    region_array = [];

    for (idx,r) in enumerate(regions)
        push!(region_array,Region(r,N,floor.(Int,load_year[idx,:])))
    end

    @info "Processing pseudoregions for VSC lines, if applicable..."
    # Load the capacity converter data, since for VSC we have to create pseudo-regions from this data
    converter_data = get_converter_capacity_data(ReEDS_data);
    #check if this dataframe is empty, so for a non-VSC case we can ignore
    VSC_append_str = "_VSC"; #necessary for naming and distinguishing VSC overlay and associated pseudoregions
    pseudoregions = [string(i)*VSC_append_str for i in converter_data[!,"r"]];
    if length(pseudoregions) > 0
        VSC_regions = [];
        for r in regions
            push!(region_array,Region(r*VSC_append_str,N,zeros(Int,N)))#should be an array of 0s for the load for pseudoregions
            push!(VSC_regions,r*VSC_append_str)#the region list must also be expanded
        end
        for r in VSC_regions
            push!(regions,r)
        end
    end
    return (VSC_append_str,N,region_array,PRAS.Regions{N,PRAS.MW}(get_name.(region_array),reduce(vcat,(get_load.(region_array)))))
end

function process_lines(ReEDS_data::CEMdata,regions::Vector{<:AbstractString},year::Int,N::Int,VSC_append_str::String)
    @info "Processing lines..."

    line_base_cap_data = get_line_capacity_data(ReEDS_data);
    ####################################################### 
    # PRAS lines
    ####################################################### 
    # Adding up line capacities between same PCA's
    @info "Processing Lines in the mapping file..."
    # Figuring out which lines belong in the PRAS System and fix the interface_regions_to, interface_regions_to indices PRAS expects

    system_line_idx = []
    region_idxs = Dict(regions .=> range(1,length(regions)))
    for (idx,row) in enumerate(eachrow(line_base_cap_data))
        pca_from = row.r
        pca_to = row.rr
        from_idx = region_idxs[pca_from]#findfirst(x->x==pca_from,regions)
        to_idx = region_idxs[pca_to]#findfirst(x->x==pca_to,regions)
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

    system_line_naming_data = DataFrames.combine(DataFrames.groupby(system_line_naming_data, ["r","rr","trtype"]), :MW => sum) #split-apply-combine b/c some lines have same name convention

    lines_array = Line[];
    for row in eachrow(system_line_naming_data)
        forward_cap = sum(line_base_cap_data[(line_base_cap_data.r.==row.r) .& (line_base_cap_data.rr.==row.rr) .& (line_base_cap_data.trtype.==row.trtype),"MW"])
        backward_cap = sum(line_base_cap_data[(line_base_cap_data.r.==row.rr) .& (line_base_cap_data.rr.==row.r) .& (line_base_cap_data.trtype.==row.trtype),"MW"])
        
        if occursin("VSC",row.trtype)
            rf = row.r*VSC_append_str #line is between pseudoregions for VSC
            rt = rt*VSC_append_str #line is between pseudoregions for VSC 
        else
            rf = row.r
            rt = row.rr
        end
        name = "$(rf)_$(rt)_$(row.trtype)"
        
        # @info "a line $name, with $forward_cap MW fwd and $backward_cap bckwrd in $(row.trtype)"
        push!(lines_array,Line(name,N,row.trtype,rf,rt,forward_cap,backward_cap,"Existing",0.0,24))#for and mttr will be defaults
    end

    #also create lines for the pseudoregions - region tx capacity
    line_base_cap_data = get_converter_capacity_data(ReEDS_data)
    for row in eachrow(line_base_cap_data) #if empty, this shouldn't iterate
        category = "VSC DC-AC converter" #this has to match the available type name in types.jl
        name = "$(row.r*VSC_append_str)_$(row.r)_$(category)"
        @info "a pseudoregion to region line $name, with $(row.MW) MW fwd and $(row.MW) bckwrd in $category"
        push!(lines_array,Line(name,N,category,row.r*VSC_append_str,row.r,row.MW,row.MW,"Existing",0.0,24))
    end

    sorted_regional_lines,temp_regions_tuple,interface_line_idxs = get_sorted_lines(lines_array,regions);

    line_forward_capacity_array = reduce(vcat,get_forward_capacity.(sorted_regional_lines));
    line_backward_capacity_array = reduce(vcat,get_backward_capacity.(sorted_regional_lines));
    λ_lines = reduce(vcat,get_λ.(sorted_regional_lines))
    μ_lines = reduce(vcat,get_μ.(sorted_regional_lines))
    new_lines = PRAS.Lines{N,1,PRAS.Hour,PRAS.MW}(get_name.(sorted_regional_lines), get_category.(sorted_regional_lines), line_forward_capacity_array, line_backward_capacity_array, λ_lines, μ_lines);

    @info "Making PRAS Interfaces..."
    new_interfaces = make_pras_interfaces(sorted_regional_lines,temp_regions_tuple,interface_line_idxs,regions);
    
    return (new_lines,new_interfaces,interface_line_idxs)
end

function split_generator_types(ReEDS_data::CEMdata,year::Int64)

    tech_types_data = get_technology_types(ReEDS_data)
    capacity_data = get_ICAP_data(ReEDS_data)

    vg_types = DataFrames.dropmissing(tech_types_data,:VRE)[:,"Column1"]
    vg_types = vg_types[vg_types .!= "csp-ns"]#csp-ns causes problems, so delete for now
    storage_types = DataFrames.dropmissing(tech_types_data,:STORAGE)[:,"Column1"]

    #clean vg/storage capacity on a regex, though there might be a better way...    
    clean_names!(vg_types)
    clean_names!(storage_types)
    vg_types = expand_types!(vg_types,15) #expand so names will match

    vg_capacity = capacity_data[(findall(in(vg_types),capacity_data.i)),:]
    storage_capacity = capacity_data[(findall(in(storage_types),capacity_data.i)),:]
    thermal_capacity = capacity_data[(findall(!in(vcat(vg_types,storage_types)),capacity_data.i)),:]
    return(thermal_capacity,storage_capacity,vg_capacity)
end

function process_thermals_with_disaggregation(thermal_builds::DataFrames.DataFrame,FOR_data::DataFrames.DataFrame,N::Int,year::Int,NEMS_path::String)#FOR_data::DataFrames.DataFrame,
    thermal_builds = thermal_builds[(thermal_builds.i.!= "csp-ns"), :] #csp-ns is not a thermal; just drop in for rn
    thermal_builds = DataFrames.combine(DataFrames.groupby(thermal_builds, ["i","r"]), :MW => sum) #split-apply-combine to handle differently vintaged entries
    EIA_db = Load_EIA_NEMS_DB(NEMS_path)
    
    all_generators = Generator[];
    native_FOR_data = FOR_data[!,"Column1"];
    lowercase_FOR_data = [lowercase(i) for i in FOR_data[!,"Column1"]];
    for row in eachrow(thermal_builds) #this loop gets the FOR for each build/tech
        #TODO: eventually, it'd be nice to lookup/pass the FOR (as dict?) and N
        @info "$(row.i) $(row.r) $(row.MW_sum) MW to disaggregate..."
        if row.i in native_FOR_data
            for_idx = findfirst(x->x==row.i,native_FOR_data)#get the idx
            gen_for = FOR_data[for_idx,"Column2"]#pull the FOR
        elseif row.i in lowercase_FOR_data
            for_idx = findfirst(x->x==row.i,lowercase_FOR_data)#get the idx
            gen_for = FOR_data[for_idx,"Column2"]#pull the FOR
        else
            gen_for = 0.05
            @info "for $(row.i) $(row.r), no gen_for is found in data, so $gen_for is used"
        end
        generator_array = disagg_existing_capacity(EIA_db,floor(Int,row.MW_sum),string(row.i),string(row.r),gen_for,N,year);
        append!(all_generators,generator_array);
    end
    return all_generators
end

function process_vg(generators_array::Vector{<:ReEDS2PRAS.Generator},vg_builds::DataFrames.DataFrame,FOR_data::DataFrames.DataFrame,ReEDS_data::CEMdata,year::Int,WeatherYear::Int,N::Int)
    #data loads
    region_mapper_df = get_region_mapping(ReEDS_data);
    cf_info = get_vg_cf_data(ReEDS_data);#load is now picked up from augur
    
    vg_profiles = cf_info["block0_values"];
    start_idx = (WeatherYear-2007)*N;

    vg_builds = DataFrames.combine(DataFrames.groupby(vg_builds, ["i","r"]), :MW => sum); #split-apply-combine to handle differently vintaged entries

    for row in eachrow(vg_builds)
        category = string(row.i)
        name = "$(category)_$(string(row.r))";
        slicer = findfirst(isequal(string(row.r)), region_mapper_df[:,"rs"]);

        if !isnothing(slicer)
            region = region_mapper_df[slicer,"*r"];
        else
            region = string(row.r)
        end
        profile_index = findfirst.(isequal.(name), (cf_info["axis0"],))[1];
        size_comp = size(vg_profiles)[1];
        profile = vg_profiles[profile_index,(start_idx+1):(start_idx+N)];
        if category in FOR_data[!,"Column1"]
            for_idx = findfirst(x->x==category,FOR_data[!,"Column1"])#get the index
            gen_for = FOR_data[for_idx,"Column2"]#pull the FOR
        else
            gen_for = .05;
        end

        if !isnothing(slicer)
            name = "$(name)_$(string(region_mapper_df[slicer,"*r"]))_" #add the region to the name
        else
            name = "$(name)_" #append for matching reasons
        end
        
        push!(generators_array,VG_Gen(name,N,region,maximum(profile),profile,category,"New",gen_for,24)) 
    end
    return generators_array
end

function process_storages(storage_builds::DataFrames.DataFrame,FOR_data::DataFrames.DataFrame,ReEDS_data::CEMdata,N::Int,regions::Vector{<:AbstractString},year::Int64)
    
    @info "handling energy capacity of storages"
    storage_energy_capacity_data = get_storage_energy_capacity_data(ReEDS_data)
    energy_capacity_df = DataFrames.combine(DataFrames.groupby(storage_energy_capacity_data, ["i","r"]), :MWh => sum) #split-apply-combine to handle differently vintaged entries

    storages_array = Storage[];
    for (idx,row) in enumerate(eachrow(storage_builds))
        name = "$(string(row.i))_$(string(row.r))"
        if string(row.i) in FOR_data[!,"Column1"]
            for_idx = findfirst(x->x==string(row.i),FOR_data[!,"Column1"])#get the idx
            gen_for = FOR_data[for_idx,"Column2"]#pull the FOR
        else
            gen_for = 0.0;
            @info "did not find FOR for storage $name $(row.r) $(row.i), so setting FOR to default value $gen_for"
        end 
        name = "$(name)_"#append for later matching
        push!(storages_array,Battery(name,N,string(row.r),string(row.i),row.MW,row.MW,energy_capacity_df[idx,"MWh_sum"],"New",1,1,1,gen_for,24)) 
    end
    storages_array, region_stor_idxs = get_sorted_components(storages_array,regions);

    stor_charge_cap_array = reduce(vcat,get_charge_capacity.(storages_array))
    stor_discharge_cap_array = reduce(vcat,get_discharge_capacity.(storages_array))
    stor_energy_cap_array = reduce(vcat,get_energy_capacity.(storages_array))
    stor_chrg_eff_array = reduce(vcat,get_charge_efficiency.(storages_array))
    stor_dischrg_eff_array = reduce(vcat,get_discharge_efficiency.(storages_array))
    stor_cryovr_eff = reduce(vcat,get_carryover_efficiency.(storages_array))
    λ_stor = reduce(vcat,get_λ.(storages_array))
    μ_stor = reduce(vcat,get_μ.(storages_array))
    return (region_stor_idxs,PRAS.Storages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(get_name.(storages_array),get_type.(storages_array),
                                                stor_charge_cap_array,stor_discharge_cap_array,stor_energy_cap_array,
                                                stor_chrg_eff_array,stor_dischrg_eff_array, stor_cryovr_eff,
                                                λ_stor,μ_stor))
end

function process_genstors(regions::Vector{<:AbstractString},N::Int)
    
    gen_stors = Gen_Storage[Gen_Storage("blank_1", N, regions[1], "blank_genstor", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)] #empty for now
    sorted_gen_stors, area_genstor_idxs  = get_sorted_components(gen_stors,regions);

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
    return (new_gen_stors,area_genstor_idxs)
end