#####################################################
# Luke
# NREL
# November 2022
# ReEDS2PRAS for NTP
# Make a PRAS System from ReEDS H5s and CSVs

function process_lines(ReEDS_data::CEMdata,regions::Vector,Year::Int,N::Int)
    @info "Processing lines..."

    line_base_cap_data = get_line_capacity_data(ReEDS_data);
    ####################################################### 
    # PRAS lines
    # ToDO: HVDC converter capacity line handling
    ####################################################### 
    # Adding up line capacities between same PCA's
    @info "Processing Lines in the mapping file..."
    
    # Figuring out which lines belong in the PRAS System and fix the interface_regions_to, interface_regions_to
    # indices PRAS expects
    system_line_idx = []
    for (idx,pca_from,pca_to) in zip(range(1,length= DataFrames.nrow(line_base_cap_data)),line_base_cap_data[:,"r"],line_base_cap_data[:,"rr"])
        from_idx = findfirst(x->x==pca_from,regions)
        to_idx = findfirst(x->x==pca_to,regions)
        if (~(isnothing(from_idx)) && ~(isnothing(to_idx)))
            if (from_idx < to_idx)
                push!(system_line_idx,idx)
                # line_base_cap_data[idx,"r"] = pca_to
                # line_base_cap_data[idx,"rr"] = pca_from
            end
        end
    end
    #order is assumed preserved in splitting these dfs for now but should likely be checked
    system_line_naming_data = line_base_cap_data[system_line_idx,:]; # all line capacities

    gdf = DataFrames.groupby(system_line_naming_data, ["r","rr","trtype"]); #split-apply-combine b/c some lines have same name convention
    system_line_naming_data = DataFrames.combine(gdf, :MW => sum); 

    lines_array = [];
    # for name,category,cap,rf,rt in zip(line_names,line_categories,system_line_naming_data[!,"MW_sum"],system_line_naming_data[!,"r"],system_line_naming_data[!,"rr"])
    for idx in range(1,DataFrames.nrow(system_line_naming_data))
        
        category = system_line_naming_data[idx,"trtype"];
        rf = system_line_naming_data[idx,"r"];
        rt = system_line_naming_data[idx,"rr"];
        # name = system_line_naming_data[idx,"r"].*"_".*system_line_naming_data[idx,"rr"].*"_".*system_line_naming_data[idx,"trtype"];
        name = rf*"_"*rt*"_"*category;
        
        forward_cap = sum(line_base_cap_data[(line_base_cap_data.r.==rf) .& (line_base_cap_data.rr.==rt) .& (line_base_cap_data.trtype.==category),"MW"]);#system_line_naming_data[idx,"MW_sum"];
        backward_cap = sum(line_base_cap_data[(line_base_cap_data.r.==rt) .& (line_base_cap_data.rr.==rf) .& (line_base_cap_data.trtype.==category),"MW"]);
        @info "a line $name, $idx with $forward_cap MW fwd and $backward_cap bckwrd"
        if forward_cap==0.0
            forward_cap=1; #capacity in each dir must be nonzero
        elseif backward_cap==0.0
            backward_cap=1; #capacity in each dir must be nonzero
        end
        push!(lines_array,line(name,N,category,rf,rt,forward_cap,backward_cap,"Existing",0.0,24))#for and mttr will be defaults
    end

    #######################################################
    # PRAS Interfaces
    #######################################################
    @info "Processing Interfaces in the mapping file..."
    num_interfaces = DataFrames.nrow(system_line_naming_data);
    @info "the number of interfaces is $num_interfaces..."

    interface_line_idxs = Array{UnitRange{Int64},1}(undef,num_interfaces);
    start_id = Array{Int64}(undef,num_interfaces); 
    for i in 1: num_interfaces
        i==1 ? start_id[i] = 1 : start_id[i] =start_id[i-1]+1
        interface_line_idxs[i] = range(start_id[i], length=1)
    end

    interface_regions_from = [findfirst(x->x==system_line_naming_data[i,"r"],regions) for i in 1:num_interfaces];
    interface_regions_to = [findfirst(x->x==system_line_naming_data[i,"rr"],regions) for i in 1:num_interfaces];
    #######################################################
    # Collecting all information to make PRAS Interfaces
    #######################################################
    @info "Making PRAS Interfaces..."
    interface_forward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
    interface_backward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
    interface_forward_capacity_array = reduce(vcat,get_forward_capacity.(lines_array));
    interface_backward_capacity_array = reduce(vcat,get_backward_capacity.(lines_array));
    # interface_forward_capacity_array = mapreduce(permutedims, vcat, get_forward_capacity.(lines_array));
    # interface_backward_capacity_array = mapreduce(permutedims, vcat, get_backward_capacity.(lines_array));
    sizeval  = size(interface_backward_capacity_array)
    @info "size of fwd array is $sizeval"
    new_interfaces = PRAS.Interfaces{N,PRAS.MW}(interface_regions_from, interface_regions_to, interface_forward_capacity_array, interface_backward_capacity_array);

    return (lines_array,new_interfaces,interface_line_idxs)
end

function split_generator_types(ReEDS_data::CEMdata,Year::Int64)

    tech_types_data = get_technology_types(ReEDS_data);
    capacity_data = get_ICAP_data(ReEDS_data);#DataFrames.DataFrame(CSV.File(capacities));

    vg_types = tech_types_data[findall(!ismissing, tech_types_data[:,"VRE"]),"Column1"]
    deleteat!(vg_types, findall(x->x=="csp-ns",vg_types)) #csp-ns causes problems, so delete for now

    storage_types = tech_types_data[findall(!ismissing, tech_types_data[:,"STORAGE"]),"Column1"]

    #clean vg/storage capacity on a regex, though there might be a better way...    
    vg_types = clean_names(vg_types)
    storage_types = clean_names(storage_types)
    vg_types = expand_types(vg_types,15) #expand so names will match

    vg_capacity = capacity_data[(findall(in(vg_types),capacity_data.i)),:]
    storage_capacity = capacity_data[(findall(in(storage_types),capacity_data.i)),:]
    thermal_capacity = capacity_data[(findall(!in(vcat(vg_types,storage_types)),capacity_data.i)),:]

    return(thermal_capacity,storage_capacity,vg_capacity)
end

function process_thermals_with_disaggregation(thermal_builds::DataFrames.DataFrame,FOR_data::DataFrames.DataFrame,N::Int,Year::Int,NEMS_path::String)#FOR_data::DataFrames.DataFrame,
    thermal_builds = thermal_builds[(thermal_builds.i.!= "csp-ns"), :] #csp-ns is not a thermal; just drop in for rn
    gdf = DataFrames.groupby(thermal_builds, ["i","r"]) #split-apply-combine to handle differently vintaged entries
    thermal_builds = DataFrames.combine(gdf, :MW => sum);
    EIA_db = Load_EIA_NEMS_DB(NEMS_path) #for now, though this is bad practice
    all_generators = [];
    native_FOR_data = FOR_data[!,"Column1"];
    lowercase_FOR_data = [lowercase(i) for i in FOR_data[!,"Column1"]];
    for (i,r,MW) in zip(thermal_builds[!,"i"],thermal_builds[!,"r"],thermal_builds[!,"MW_sum"])
        #eventually, it'd be nice to lookup/pass the FOR and N
        @info "$i $r $MW MW to disaggregate..."
        #get the FOR
        if i in native_FOR_data
            for_idx = findfirst(x->x==i,native_FOR_data)#get the idx
            gen_for = FOR_data[for_idx,"Column2"]#pull the FOR
        elseif i in lowercase_FOR_data
            for_idx = findfirst(x->x==i,lowercase_FOR_data)#get the idx
            gen_for = FOR_data[for_idx,"Column2"]#pull the FOR
        else
            gen_for = 0.05
            @info "for $i $r, no gen_for is found in data, so $gen_for is used"
            # error("we really should always find a gen_for, but did not for $i $r") #default val
        end
        generator_array = disagg_existing_capacity(EIA_db,floor.(Int,MW),string.(i),string.(r),gen_for,N,Year);
        append!(all_generators,generator_array);
    end
    return all_generators
end

function process_vg(generators_array::Vector,vg_builds::DataFrames.DataFrame,FOR_data::DataFrames.DataFrame,ReEDS_data::CEMdata,Year::Int,WeatherYear::Int,N::Int)
    #data loads
    # region_mapper = joinpath(ReEDSfilepath,"inputs_case","rsmap.csv");
    region_mapper_df = get_region_mapping(ReEDS_data)#DataFrames.DataFrame(CSV.File(region_mapper));
    cf_info = get_vg_cf_data(ReEDS_data);#; #load is now picked up from augur
    
    vg_profiles = cf_info["block0_values"];
    start_idx = (WeatherYear-2007)*N;

    #split-apply-combine to group...
    gdf = DataFrames.groupby(vg_builds, ["i","r"]); #split-apply-combine to handle differently vintaged entries
    vg_builds = DataFrames.combine(gdf, :MW => sum);
    # @info "vg builds are now $vg_builds"

    for idx in range(1,DataFrames.nrow(vg_builds))
        category = string(vg_builds[idx,"i"]);
        name = category*"_"*string(vg_builds[idx,"r"]);
        slicer = findfirst(isequal(string(vg_builds[idx,"r"])), region_mapper_df[:,"rs"]);
        if !isnothing(slicer)
            region = region_mapper_df[slicer,"*r"];
        else
            region = string(vg_builds[idx,"r"])
        end
        profile_index = findfirst.(isequal.(name), (cf_info["axis0"],))[1];
        size_comp = size(vg_profiles)[1];
        profile = vg_profiles[profile_index,(start_idx+1):(start_idx+N)];
        if category in FOR_data[!,"Column1"]
            for_idx = findfirst(x->x==category,FOR_data[!,"Column1"])#get the idx
            gen_for = FOR_data[for_idx,"Column2"]#pull the FOR
        else
            gen_for = .05;
        end

        if !isnothing(slicer)
            name = name*"_"*string(region_mapper_df[slicer,"*r"])*"_"; #add the region to the name
        else
            name = name*"_"; #append for matching reasons
        end
        
        push!(generators_array,vg_gen(name,N,region,maximum(profile),profile,category,"New",gen_for,24)) 
    end
    return generators_array
end

function process_storages(storage_builds::DataFrames.DataFrame,FOR_data::DataFrames.DataFrame,ReEDS_data::CEMdata,N::Int,regions::Vector,Year::Int64)
    #get the vector of appropriate indices
    @info "handling power capacity of storages"
    storage_names = [string(storage_builds[!,"i"][i])*"_"*string(storage_builds[!,"r"][i]) for i=1:DataFrames.nrow(storage_builds)];
    storage_capacities = storage_builds[!,"MW"];#want capacities
    #get regions and convert them to appropriate data where relevant
    storage_categories = string.(storage_builds[!,"i"]);#[string(vg_builds[!,"i"][i]) for i=1:DataFrames.nrow(vg_builds)]
    storage_regions = string.(storage_builds[!,"r"]);
    region_stor_capacity_idxs,stor_capacity_region_order = sort_gens(storage_regions,regions,length(regions))
    
    @info "handling energy capacity of storages"
    storage_energy_capacity_data = get_storage_energy_capacity_data(ReEDS_data)
       
    gdf = DataFrames.groupby(storage_energy_capacity_data, ["i","r"]) #split-apply-combine to handle differently vintaged entries
    annual_data_sum = DataFrames.combine(gdf, :MWh => sum);
    storage_energy_capacity_regions = string.(annual_data_sum[!,"r"]);
    #reorder to ensure proper order!
    region_stor_idxs,stor_region_order = sort_gens(storage_energy_capacity_regions,regions,length(regions))
    
    #To-do, this can almost certainly be refactored
    storages_array = [];
    for (name,region,category,capacity,energy_capacity) in zip(storage_names[stor_capacity_region_order],storage_regions,storage_categories[stor_capacity_region_order],storage_capacities[stor_capacity_region_order],annual_data_sum[stor_region_order,"MWh_sum"])
        if category in FOR_data[!,"Column1"]
            for_idx = findfirst(x->x==category,FOR_data[!,"Column1"])#get the idx
            gen_for = FOR_data[for_idx,"Column2"]#pull the FOR
        else
            gen_for = 0.0;
            @info "did not find FOR for storage $name $region $category, so setting FOR to $gen_for"
        end 
        name = name*"_";#append for matching
        push!(storages_array,battery(name,N,region,category,capacity,capacity,energy_capacity,"New",1,1,1,gen_for,24)) 
    end
    return (storages_array,region_stor_idxs)
end

function make_pras_system_from_mapping_info(ReEDSfilepath::String, Year::Int64, CaseLabel::String, NEMS_path::String)
    #######################################################
    # Loading the necessary mapping files and data
    #######################################################
    #check validity of input weather and ReEDS year
    WEATHERYEAR = 2012 #just for now, force this
    if !isfile(joinpath(ReEDSfilepath,"ReEDS_Augur","augur_data","plot_load_"*string(Year)*".h5"))
        error("The year $Year does not have an associated Augur file (only checking load file for now). Are you sure ReeDS was run and Augur results saved for $Year?")
    end
    if WEATHERYEAR ∉ [2007,2008,2009,2010,2011,2012,2013] 
        error("The weather year $WEATHERYEAR is not a valid VG profile year. Should be an Int in 2007-2013 currently")
    end
    # Load the mapping metadata JSON file
    @info "Fetching ReEDS case data to build PRAS System..."
    ReEDS_data = ReEDSdata(ReEDSfilepath,Year);
    load_info = get_load_file(ReEDS_data);#HDF5.h5read(joinpath(ReEDSfilepath,"ReEDS_Augur","augur_data","plot_load_"*string(Year)*".h5"),"data"); #load is now picked up from augur
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
        push!(region_array,region(r,N,floor.(Int,load_year[idx,:])))
    end

    new_regions = PRAS.Regions{N,PRAS.MW}(get_name.(region_array),reduce(vcat,(get_load.(region_array))));

    #######################################################
    # PRAS Region Gen Index 
    # **TODO: Should 0 MW generators be allowed after disaggregation?
    # **TODO: Should hydro be split out as a generator-storage?
    # **TODO: is it important to also handle planned outages?
    #######################################################
    line_array,new_interfaces,interface_line_idxs = process_lines(ReEDS_data,regions,Year,8760);

    line_forward_capacity_array = reduce(vcat,get_forward_capacity.(line_array));#permutedims(hcat(get_forward_capacity.(line_array)...));
    line_backward_capacity_array = reduce(vcat,get_backward_capacity.(line_array));#permutedims(hcat(get_backward_capacity.(line_array)...));
    λ_lines = reduce(vcat,get_λ.(line_array))#permutedims(hcat(get_λ.(line_array)...));
    μ_lines = reduce(vcat,get_μ.(line_array))#permutedims(hcat(get_μ.(line_array)...));
    new_lines = PRAS.Lines{N,1,PRAS.Hour,PRAS.MW}(get_name.(line_array), get_category.(line_array), line_forward_capacity_array, line_backward_capacity_array, λ_lines, μ_lines);
    
    @info "splitting thermal, storage, vg generator types from installed ReEDS capacities..."
    thermal,storage,vg = split_generator_types(ReEDS_data,Year);

    @info "reading in ReEDS generator-type forced outage data..."
    forced_outage_data = get_forced_outage_data(ReEDS_data);

    @info "reading thermals..."
    gens = process_thermals_with_disaggregation(thermal,forced_outage_data,N,Year,NEMS_path);
    @info "reading vg..."
    gens = process_vg(gens,vg,forced_outage_data,ReEDS_data,Year,WEATHERYEAR,N);
    @info "...vg read, translating to PRAS gens"
    area_gen_idxs,region_order = sort_gens(get_region.(gens),regions,length(regions))
    gens = gens[region_order]; #now they are properly sorted, we hope!
    
    capacity_matrix = reduce(vcat,get_capacity.(gens));
    λ_matrix = reduce(vcat,get_λ.(gens));
    μ_matrix = reduce(vcat,get_μ.(gens));

    new_generators = PRAS.Generators{N,1,PRAS.Hour,PRAS.MW}(get_name.(gens),get_type.(gens),capacity_matrix,λ_matrix,μ_matrix);

    #######################################################
    # PRAS Timestamps
    #######################################################
    @info "Processing PRAS timestamps..."
    start_datetime = Dates.DateTime(string(Year)*"-01-01");
    finish_datetime = start_datetime + Dates.Hour(N-1);
    my_timestamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);

    #######################################################
    # PRAS Storages
    #######################################################
    @info "Processing Storages..."
    storages,area_stor_idxs = process_storages(storage,forced_outage_data,ReEDS_data,N,regions,Year)

    stor_charge_cap_array = reduce(vcat,get_charge_capacity.(storages));#permutedims(hcat(get_charge_capacity.(storages)...))
    stor_discharge_cap_array = reduce(vcat,get_discharge_capacity.(storages));
    stor_energy_cap_array = reduce(vcat,get_energy_capacity.(storages));
    stor_chrg_eff_array = reduce(vcat,get_charge_efficiency.(storages));
    stor_dischrg_eff_array = reduce(vcat,get_discharge_efficiency.(storages));
    stor_cryovr_eff = reduce(vcat,get_carryover_efficiency.(storages));
    λ_stor = reduce(vcat,get_λ.(storages));
    μ_stor = reduce(vcat,get_μ.(storages));
    new_storage = PRAS.Storages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(get_name.(storages),get_type.(storages),
                                                stor_charge_cap_array,stor_discharge_cap_array,stor_energy_cap_array,
                                                stor_chrg_eff_array,stor_dischrg_eff_array, stor_cryovr_eff,
                                                λ_stor,μ_stor);
    #######################################################
    # PRAS GeneratorStorages
    # No GeneratorStorages in this system
    #######################################################
    @info "Processing GeneratorStorages [EMPTY FOR NOW].."
    gen_stor_names = String[];
    gen_stor_categories = String[];

    n_genstors = 0;

    gen_stor_charge_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_discharge_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_enrgy_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_charge_eff = Matrix{Float64}(undef, n_genstors, N);
    gen_stor_discharge_eff = Matrix{Float64}(undef, n_genstors, N);
    gen_stor_cryovr_eff = Matrix{Float64}(undef, n_genstors, N);
    gen_stor_inflow_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_gridwdr_cap_array = Matrix{Int64}(undef, n_genstors, N);
    gen_stor_gridinj_cap_array = Matrix{Int64}(undef, n_genstors, N);

    λ_genstors = Matrix{Float64}(undef, n_genstors, N);
    μ_genstors = Matrix{Float64}(undef, n_genstors, N);

    new_gen_stors = PRAS.GeneratorStorages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(gen_stor_names,gen_stor_categories,
                                                            gen_stor_charge_cap_array, gen_stor_discharge_cap_array, gen_stor_enrgy_cap_array,
                                                            gen_stor_charge_eff, gen_stor_discharge_eff, gen_stor_cryovr_eff,
                                                            gen_stor_inflow_array, gen_stor_gridwdr_cap_array, gen_stor_gridinj_cap_array,
                                                            λ_genstors, μ_genstors);
                                                            
    area_genstor_idxs = fill(1:0, num_areas);
    pras_system = PRAS.SystemModel(new_regions, new_interfaces, new_generators, area_gen_idxs, new_storage, area_stor_idxs, new_gen_stors,
                            area_genstor_idxs, new_lines,interface_line_idxs,my_timestamps);
    #save PRAS system somewhere we can use it?
    PRAS.savemodel(pras_system,joinpath(ReEDSfilepath,"outputs",CaseLabel*".pras"))
    short,flow = run_pras_system(pras_system,2)#just two for now to save time but eventually more
    return(pras_system)
    
end