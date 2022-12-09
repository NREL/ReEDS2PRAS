#####################################################
# Luke
# NREL
# November 2022
# ReEDS2PRAS for NTP
# Make a PRAS System from ReEDS H5s and CSVs

function process_lines(ReEDSfilepath::String,regions::Vector,Year::Int,N::Int)
    @info "Processing lines..."

    if isfile(joinpath(ReEDSfilepath,"outputs","tran_out.csv"))
        line_base_cap_csv_location = joinpath(ReEDSfilepath,"outputs","tran_out.csv")
        line_base_cap_data = DataFrames.DataFrame(CSV.File(line_base_cap_csv_location));#load the csv
        DataFrames.rename!(line_base_cap_data, (names(line_base_cap_data) .=> ["*r","rr","trtype","Year","MW"])...)#rename the columns
        line_base_cap_data = line_base_cap_data[line_base_cap_data.Year.==Year,:] #filter by year, ReEDS line capacity is cumulative 

    elseif isfile(joinpath(ReEDSfilepath,"inputs_case","trancap_init_energy.csv"))
        @info "---> on error, read old file structure for now, but, generally, this is a cause for concern..."
        line_cap_base_csv_location = joinpath(ReEDSfilepath,"inputs_case","trancap_init_energy.csv") #there is also a _prm file, not sure which is right?
        line_cap_additional_csv_location = joinpath(ReEDSfilepath,"inputs_case","trancap_fut.csv")

        line_base_cap_data = DataFrames.DataFrame(CSV.File(line_cap_base_csv_location;select=["*r","rr","trtype","MW"]));
        line_additional_cap_data = DataFrames.DataFrame(CSV.File(line_cap_additional_csv_location;select=["*r","rr","trtype","MW"]));
    else
        @info "---> default file read for transmission capacity fails, try reading alternative file structure"
        line_cap_base_csv_location = joinpath(ReEDSfilepath,"inputs_case","trancap_init.csv"); #there is also a _prm file, not sure which is right?
        line_base_cap_data = DataFrames.DataFrame(CSV.File(line_cap_base_csv_location;select=["*r","rr","trtype","MW"]));
    end
    ####################################################### 
    # PRAS lines
    # ToDO: Do we need to handle HVDC lines differently?
    # ToDO: file naming convention on ReEDS side is tough to handle - strings must be exactly right. Prefer to pass that in
    ####################################################### 
    # Adding up line capacities between same PCA's
    @info "Processing Lines in the mapping file..."
    # base_cap_pca_pairs = line_base_cap_data[!,"*r"].*"_".*line_base_cap_data[!,"rr"];
    # unique_base_cap_pca_pairs  = unique(base_cap_pca_pairs);
    # println(base_cap_pca_pairs)
    # for pca_pair in unique_base_cap_pca_pairs
    #     pca_idxs = findall(x -> x == pca_pair, base_cap_pca_pairs);
    #     println(pca_idxs)
    #     if (length(pca_idxs) >1)
    #         sum_cap = 0
    #         for pca_idx in pca_idxs
    #             println(pca_idx)
    #             sum_cap+=line_base_cap_data[pca_idx,"MW"]
    #         end
    #         println("add")
    #         line_base_cap_data[pca_idxs[1],"MW"] = sum_cap
    #         println(line_base_cap_data)
    #         for del_idx in range(2,length=length(pca_idxs)-1)
    #             println("a del")
    #             delete!(line_base_cap_data,pca_idxs[del_idx])
    #         end
    #         println("...")
    #         println(line_base_cap_data)
    #     end
    # end
    line_base_cap_data[!, "direction"] = ones(size(line_base_cap_data)[1])
    
    # Figuring out which lines belong in the PRAS System and fix the interface_regions_to, interface_regions_to
    # indices PRAS expects
    system_line_idx = []
    for (idx,pca_from,pca_to) in zip(range(1,length= DataFrames.nrow(line_base_cap_data)),line_base_cap_data[:,"*r"],line_base_cap_data[:,"rr"])
        from_idx = findfirst(x->x==pca_from,regions)
        to_idx = findfirst(x->x==pca_to,regions)
        if (~(isnothing(from_idx)) && ~(isnothing(to_idx)))
            push!(system_line_idx,idx)
            if (from_idx > to_idx)
                line_base_cap_data[idx,"*r"] = pca_to
                line_base_cap_data[idx,"rr"] = pca_from
                line_base_cap_data[idx,"direction"] = 2.
            #     push!(back_idx,idx)
            # elseif (to_idx > from_idx)
            #     push!(fwd_idx,idx)
            end
        end
    end
    #order is assumed preserved in splitting these dfs for now but should likely be checked
    system_line_naming_data = line_base_cap_data[system_line_idx,:]; # all line capacities
    # println(system_line_naming_data)
    # line_idxs,line_region_order = sort_gens(system_line_naming_data[!,"*r"],regions,length(regions))
    # system_line_naming_data = system_line_naming_data[line_region_order,:]

    gdf = DataFrames.groupby(system_line_naming_data, ["*r","rr","trtype"]) #split-apply-combine b/c some lines have same name convention
    system_line_naming_data = DataFrames.combine(gdf, :MW => sum); 

    line_names = system_line_naming_data[!,"*r"].*"_".*system_line_naming_data[!,"rr"].*"_".*system_line_naming_data[!,"trtype"];
    line_categories = string.(system_line_naming_data[!,"trtype"]);

    #######################################################
    # Collecting all information to make PRAS Lines
    #######################################################
    n_lines = DataFrames.nrow(system_line_naming_data);
    line_forward_capacity_array = Matrix{Int64}(undef, n_lines, N);
    line_backward_capacity_array = Matrix{Int64}(undef, n_lines, N);

    λ_lines = Matrix{Float64}(undef, n_lines, N); 
    μ_lines = Matrix{Float64}(undef, n_lines, N); 
    for (idx,line_cap) in enumerate(system_line_naming_data[!,"MW_sum"])
        line_forward_capacity_array[idx,:] = fill.(floor.(Int,line_cap),1,N); #forward and backward capacities are the same
        line_backward_capacity_array[idx,:] = fill.(floor.(Int,line_cap),1,N);

        λ_lines[idx,:] .= 0.0; 
        μ_lines[idx,:] .= 1.0; 
    end

    new_lines = PRAS.Lines{N,1,PRAS.Hour,PRAS.MW}(line_names, line_categories, line_forward_capacity_array, line_backward_capacity_array, λ_lines, μ_lines);

    #######################################################
    # PRAS Interfaces
    #######################################################
    @info "Processing Interfaces in the mapping file..."
    num_interfaces = DataFrames.nrow(system_line_naming_data);

    interface_line_idxs = Array{UnitRange{Int64},1}(undef,num_interfaces);
    start_id = Array{Int64}(undef,num_interfaces); 
    for i in 1: num_interfaces
        i==1 ? start_id[i] = 1 : start_id[i] =start_id[i-1]+1
        interface_line_idxs[i] = range(start_id[i], length=1)
    end

    interface_regions_from = [findfirst(x->x==system_line_naming_data[i,"*r"],regions) for i in 1:num_interfaces];
    interface_regions_to = [findfirst(x->x==system_line_naming_data[i,"rr"],regions) for i in 1:num_interfaces];
    #######################################################
    # Collecting all information to make PRAS Interfaces
    #######################################################
    @info "Making PRAS Interfaces..."
    interface_forward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
    interface_backward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
    for i in 1:num_interfaces
        interface_forward_capacity_array[i,:] .=  line_forward_capacity_array[i,:]
        interface_backward_capacity_array[i,:] =  line_backward_capacity_array[i,:]
    end

    new_interfaces = PRAS.Interfaces{N,PRAS.MW}(interface_regions_from, interface_regions_to, interface_forward_capacity_array, interface_backward_capacity_array);

    return (new_lines,new_interfaces,interface_line_idxs)
end

function split_generator_types(ReEDSfilepath::String,Year::Int64)
    tech_types = joinpath(ReEDSfilepath,"inputs_case","tech-subset-table.csv");
    tech_types_data = DataFrames.DataFrame(CSV.File(tech_types));

    vg_types = tech_types_data[findall(!ismissing, tech_types_data[:,"VRE"]),"Column1"]
    deleteat!(vg_types, findall(x->x=="csp-ns",vg_types)) #csp-ns causes problems, so delete for now

    storage_types = tech_types_data[findall(!ismissing, tech_types_data[:,"STORAGE"]),"Column1"]

    #clean vg/storage capacity on a regex, though there might be a better way...    
    vg_types = clean_names(vg_types)
    storage_types = clean_names(storage_types)
    vg_types = expand_types(vg_types,15) #expand so names will match

    capacities = joinpath(ReEDSfilepath,"outputs","cap.csv")
    capacity_data = DataFrames.DataFrame(CSV.File(capacities))
    annual_data = capacity_data[capacity_data.Dim3.==Year,:] #only built capacity in the PRAS-year, since ReEDS capacity is cumulative    

    vg_capacity = annual_data[(findall(in(vg_types),annual_data.Dim1)),:]
    storage_capacity = annual_data[(findall(in(storage_types),annual_data.Dim1)),:]
    thermal_capacity = annual_data[(findall(!in(vcat(vg_types,storage_types)),annual_data.Dim1)),:]

    return(thermal_capacity,storage_capacity,vg_capacity)
end

function create_generators_from_data!(gen_matrix,gen_names,gen_categories,FOR_data)
    #grab the FORS
    n_gen,N = size(gen_matrix)[1],size(gen_matrix)[2];
    gen_cap_array = Matrix{Int64}(undef, n_gen, N);
    λ_gen = Matrix{Float64}(undef, n_gen, N);
    μ_gen = Matrix{Float64}(undef, n_gen, N);
    for idx in 1:n_gen
        gen_cap_array[idx,:] = floor.(Int,gen_matrix[idx,:]); #converts to int
        
        #use conditional to set failure and recovery probabilities for generators
        #eventually this will have to call to Sinnott's data & be much more sophisticated
        if gen_categories[idx] in FOR_data[!,"Column1"]
            for_idx = findfirst(x->x==gen_categories[idx],FOR_data[!,"Column1"])#get the idx
            gen_for = FOR_data[for_idx,"Column2"]#pull the FOR
            lambda,mu = FOR_to_transitionprobs(gen_for)#run the helper fxn
            λ_gen[idx,:] .= lambda;
            μ_gen[idx,:] .= mu;
        else
            λ_gen[idx,:] .= 0.01; 
            μ_gen[idx,:] .= 0.2;
        end
    end


    return (gen_names,gen_categories,gen_cap_array,λ_gen,μ_gen);
end

function process_thermals(thermal_builds::DataFrames.DataFrame,FOR_data::DataFrames.DataFrame,N::Int)
    thermal_builds = thermal_builds[(thermal_builds.Dim1.!= "csp-ns"), :] #csp-ns is not a thermal; just drop in for rn
    #get the vector of appropriate indices
    thermal_names = [string(thermal_builds[!,"Dim1"][i])*"_"*string(thermal_builds[!,"Dim2"][i]) for i=1:DataFrames.nrow(thermal_builds)]
    thermal_capacities = thermal_builds[!,"Val"]#want capacities

    #get regions and convert them to appropriate data where relevant
    # thermal_categories = [string(thermal_builds[!,"Dim1"][i]) for i=1:DataFrames.nrow(thermal_builds)]
    thermal_categories = string.(thermal_builds[!,"Dim1"])
    thermal_regions = string.(thermal_builds[!,"Dim2"])

    thermal_cap_factors = ones(length(thermal_names),N)

    return(create_generators_from_data!(thermal_capacities.*thermal_cap_factors,thermal_names,thermal_categories,FOR_data),thermal_regions)
end

function process_vg(vg_builds::DataFrames.DataFrame,FOR_data::DataFrames.DataFrame,ReEDSfilepath::String,Year::Int,WeatherYear::Int,N::Int)
    #get the vector of appropriate indices
    vg_names = [string(vg_builds[!,"Dim1"][i])*"_"*string(vg_builds[!,"Dim2"][i]) for i=1:DataFrames.nrow(vg_builds)];
    vg_capacities = vg_builds[!,"Val"];#want capacities

    #get regions and convert them to appropriate data where relevant
    vg_categories = string.(vg_builds[!,"Dim1"]);#[string(vg_builds[!,"Dim1"][i]) for i=1:DataFrames.nrow(vg_builds)]
    vg_regions = string.(vg_builds[!,"Dim2"]);
    
    #pull in mapping info for the regions, since vg is originally not assigned by region
    region_mapper = joinpath(ReEDSfilepath,"inputs_case","rsmap.csv");
    region_mapper_df = DataFrames.DataFrame(CSV.File(region_mapper));
    for (idx,region) in enumerate(vg_regions)
        slicer = findfirst(isequal(region), region_mapper_df[:,"rs"]);
        if !isnothing(slicer)
            vg_regions[idx] = region_mapper_df[slicer,"*r"]
        end
    end

    #load in vg data from augur
    # cf_info = HDF5.h5read(joinpath(ReEDSfilepath,"inputs_case","recf.h5"), "data")
    cf_info = HDF5.h5read(joinpath(ReEDSfilepath,"ReEDS_Augur","augur_data","plot_vre_gen_"*string(Year)*".h5"),"data"); #load is now picked up from augur
    indices = findfirst.(isequal.(vg_names), (cf_info["axis0"],));
    # indices = findfirst.(isequal.(vg_names), (cf_info["block0_items"],))
    vg_profiles = cf_info["block0_values"];
    start_idx = (WeatherYear-2007)*N;
    # retained_cap_factors = cap_factors[indices,1:N] #slice recf, just take first year for now until more is known
    retained_vg_profiles = vg_profiles[indices,(start_idx+1):(start_idx+N)]; #return the profiles, no longer unitized
    return(create_generators_from_data!(retained_vg_profiles,vg_names,vg_categories,FOR_data),vg_regions)
    # return(create_generators_from_data!(vg_capacities.*retained_cap_factors,vg_names,vg_categories,FOR_data),vg_regions)
end

function process_storages(storage_builds::DataFrames.DataFrame,FOR_data::DataFrames.DataFrame,ReEDSfilepath::String,N::Int,regions::Vector,Year::Int64)
    #get the vector of appropriate indices
    @info "handling power capacity of storages"
    storage_names = [string(storage_builds[!,"Dim1"][i])*"_"*string(storage_builds[!,"Dim2"][i]) for i=1:DataFrames.nrow(storage_builds)];
    storage_capacities = storage_builds[!,"Val"];#want capacities
    storage_cap_factors = ones(length(storage_names),N);
    #get regions and convert them to appropriate data where relevant
    storage_categories = string.(storage_builds[!,"Dim1"]);#[string(vg_builds[!,"Dim1"][i]) for i=1:DataFrames.nrow(vg_builds)]
    storage_regions = string.(storage_builds[!,"Dim2"]);
    region_stor_capacity_idxs,stor_capacity_region_order = sort_gens(storage_regions,regions,length(regions))
    stor_names,stor_categories,stor_discharge_cap_array,λ_stor,μ_stor= create_generators_from_data!(storage_capacities[stor_capacity_region_order].*storage_cap_factors,storage_names[stor_capacity_region_order],storage_categories[stor_capacity_region_order],FOR_data);
    stor_charge_cap_array = stor_discharge_cap_array #make this equal for now
    
    @info "handling energy capacity of storages"
    storage_energy_capacities = joinpath(ReEDSfilepath,"outputs","stor_energy_cap.csv")
    storage_energy_capacity_data = DataFrames.DataFrame(CSV.File(storage_energy_capacities))
    annual_data = storage_energy_capacity_data[storage_energy_capacity_data.Dim4.==Year,:] #only built capacity in the PRAS-year, since ReEDS capacity is cumulative    
    
    gdf = DataFrames.groupby(annual_data, ["Dim1","Dim3"]) #split-apply-combine to handle differently vintaged entries
    annual_data_sum = DataFrames.combine(gdf, :Val => sum);
    storage_energy_capacity_regions = string.(annual_data_sum[!,"Dim3"]);
    #reorder to ensure proper order!
    region_stor_idxs,stor_region_order = sort_gens(storage_energy_capacity_regions,regions,length(regions))
    stor_energy_cap_array = annual_data_sum[stor_region_order,"Val_sum"].*storage_cap_factors
    stor_energy_cap_array = floor.(Int, stor_energy_cap_array)#go to int

    #now we can do the matrix math on the sorted vector of energy capacities
    # stor_energy_cap_array = stor_discharge_cap_array #this is a terrible idea, but for now
    
    #efficiency and carryover; fill all as ones for now
    stor_chrg_eff_array = ones(length(stor_names),N);
    stor_dischrg_eff_array = ones(length(stor_names),N);
    stor_cryovr_eff = ones(length(stor_names),N);

    new_storage = PRAS.Storages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(stor_names,stor_categories,
                                                stor_charge_cap_array,stor_discharge_cap_array,stor_energy_cap_array,
                                                stor_chrg_eff_array,stor_dischrg_eff_array, stor_cryovr_eff,
                                                λ_stor,μ_stor);

    return(new_storage,region_stor_idxs)
end

function make_pras_system_from_mapping_info(ReEDSfilepath::String, Year::Int64, CaseLabel::String)
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

    #load-related information
    # load_info = HDF5.h5read(joinpath(ReEDSfilepath,"inputs_case","load.h5"), "data");
    # "StandScen_MidCase_"*string(test_year)
    load_info = HDF5.h5read(joinpath(ReEDSfilepath,"ReEDS_Augur","augur_data","plot_load_"*string(Year)*".h5"),"data"); #load is now picked up from augur
    load_data = load_info["block0_values"];
    regions = load_info["block0_items"];

    #To-Do: can get a bug/error here if ReEDS case lacks multiple load years
    slicer = findfirst(isequal(WEATHERYEAR), load_info["axis1_level1"]); #2012 weather year, Gbut, in general, this should be a user input
    indices = findall(i->(i==(slicer-1)),load_info["axis1_label1"]); #I think b/c julia indexes from 1 we need -1 here
    #probably want some kind of assert/error here if indices are empty that states the year is invalid
    
    load_year = load_data[:,indices[1:8760]]; #should be regionsX8760, which is now just enforced
    
    @info "Processing Areas in the mapping file into PRAS regions..."
    # area_names = string.(mapping_data[!,"PlainInter"]);
    num_areas = length(regions); 
    N = size(load_year)[2];#DataFrames.DataFrames.nrow(load_data)
    region_load = Array{Int64,2}(undef,num_areas,N);

    for (idx,region) in enumerate(regions)
        region_load[idx,:]=floor.(Int,load_year[idx,:]); #converts to int
    end
    new_regions = PRAS.Regions{N,PRAS.MW}(regions, region_load);

    #######################################################
    # PRAS Region Gen Index 
    # **TODO: Figure out how to not take in account "0.0" gen_mw_max generators!!
    # **TODO : If fuel info is available, figure out how to not count them to avoid double counting.
    #######################################################
    new_lines,new_interfaces,interface_line_idxs = process_lines(ReEDSfilepath,regions,Year,8760);

    #generation capacity
    @info "splitting thermal, storage, vg generator types from installed ReEDS capacities..."
    thermal,storage,vg = split_generator_types(ReEDSfilepath,Year);

    #pull the generators into appropriate STRUCTS, with defaults

    @info "reading in ReEDS generator-type forced outage data..."
    # is it important to also handle planned outages?
    forced_outage_path = joinpath(ReEDSfilepath,"inputs_case","outage_forced.csv") #there is also a _prm file, not sure which is right?
    forced_outage_data = DataFrames.DataFrame(CSV.File(forced_outage_path,header=false));

    #for thermal, we need outage stuff
    #for thermal, we will need a disaggreggation helper fxn. Possibly a couple. May need EIA860 data
    #for vg, we need profiles
    vg_tup = process_vg(vg,forced_outage_data,ReEDSfilepath,Year,WEATHERYEAR,N);
    thermal_tup = process_thermals(thermal,forced_outage_data,N);

    gen_tup = [vcat(a,b) for (a,b) in zip(thermal_tup[1],vg_tup[1])];
    area_gen_idxs,region_order = sort_gens(vcat(thermal_tup[2],vg_tup[2]),regions,length(regions))
    new_generators = PRAS.Generators{N,1,PRAS.Hour,PRAS.MW}(gen_tup[1][region_order],gen_tup[2][region_order],gen_tup[3][region_order,:],gen_tup[4][region_order,:],gen_tup[5][region_order,:]); #there's probably a way to better unpack this

    #######################################################
    # PRAS Timestamps
    #######################################################
    @info "Processing PRAS timestamps..."
    start_datetime = Dates.DateTime(string(Year)*"-01-01");
    finish_datetime = start_datetime + Dates.Hour(N-1);
    my_timestamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);

    #######################################################
    # PRAS Storages
    # Storages are added 11.28
    #######################################################
    @info "Processing Storages..."
    new_storage,area_stor_idxs = process_storages(storage,forced_outage_data,ReEDSfilepath,N,regions,Year)
    
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
    short,flow = run_pras_system(pras_system,10)
    return(pras_system)
    
end