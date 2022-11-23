
#####################################################
# Luke
# NREL
# November 2022
# ReEDS2PRAS for NTP
# Make a PRAS System from ReEDS H5s and CSVs

function process_lines(ReEDSfilepath::String)
    @info "Processing lines..."
    line_cap_base_csv_location = joinpath(ReEDSfilepath,"inputs_case","trancap_init_energy.csv") #there is also a _prm file, not sure which is right?
    line_cap_additional_csv_location = joinpath(ReEDSfilepath,"inputs_case","trancap_fut.csv")

    line_base_cap_data = DataFrames.DataFrame(CSV.File(line_cap_base_csv_location;select=["*r","rr","trtype","MW"]));
    line_additional_cap_data = DataFrames.DataFrame(CSV.File(line_cap_additional_csv_location;select=["*r","rr","trtype","MW"]));
    
    ####################################################### 
    # PRAS lines
    # ToDO: Do we need to handle HVDC lines differently?
    # ToDO: file naming convention on ReEDS side is tough to handle - strings must be exactly right. Prefer to pass that in
    ####################################################### 
    # Adding up line capacities between same PCA's
    @info "Processing Lines in the mapping file..."
    base_cap_pca_pairs = line_base_cap_data[!,"*r"].*"_".*line_base_cap_data[!,"rr"];
    unique_base_cap_pca_pairs  = unique(base_cap_pca_pairs);
    
    for pca_pair in unique_base_cap_pca_pairs
        pca_idxs = findall(x -> x == pca_pair, base_cap_pca_pairs);
        if (length(pca_idxs) >1)
            sum_cap = 0
            for pca_idx in pca_idxs
                sum_cap+=line_base_cap_data[pca_idx,"MW"]
            end
            line_base_cap_data[pca_idxs[1],"MW"] = sum_cap
            for del_idx in range(2,length=length(pca_idxs)-1)
                delete!(line_base_cap_data,pca_idxs[del_idx])
            end
        end
    end
    # Adding additional capacities to base capacities
    # for (idx,pca_from,pca_to,add_cap) in zip(range(1, length= DataFrames.nrow(line_additional_cap_data)),line_additional_cap_data[:,"r"],line_additional_cap_data[:,"rr"],line_additional_cap_data[:,"MW"])
    #     temp_idx = findall(x->x==pca_from, line_base_cap_data[!,"r"])
    #     map_idx  = findall(x->x==pca_to, line_base_cap_data[temp_idx,"rr"])
    #     if (length(map_idx) !=0)
    #         line_base_cap_data[temp_idx[map_idx][1],"MW"] = line_base_cap_data[temp_idx[map_idx][1],"MW"] + add_cap
    #     else
    #         push!(line_base_cap_data,line_additional_cap_data[idx,:])
    #     end
    # end
    return (line_base_cap_data)
end

function split_generator_types(ReEDSfilepath::String,Year::Int64)
    tech_types = joinpath(ReEDSfilepath,"inputs_case","tech-subset-table.csv")
    tech_types_data = DataFrames.DataFrame(CSV.File(tech_types));

    vg_types = tech_types_data[findall(!ismissing, tech_types_data[:,"VRE"]),"Column1"]
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


function create_generators_from_data(gen_matrix,gen_names,gen_categories,regions)
    n_gen,N = size(gen_matrix)[1],size(gen_matrix)[2];
    gen_cap_array = Matrix{Int64}(undef, n_gen, N);
    λ_gen = Matrix{Float64}(undef, n_gen, N);
    μ_gen = Matrix{Float64}(undef, n_gen, N);
    for idx in 1:n_gen
        gen_cap_array[idx,:] = floor.(Int,gen_matrix[idx,:]); #converts to int
        λ_gen[idx,:] .= 0.0; #evetually this will have to call to Sinnott's data 
        μ_gen[idx,:] .= 1.0;
    end


    return (gen_names, gen_categories, gen_cap_array , λ_gen ,μ_gen);
end

function process_thermals(thermal_builds::DataFrames.DataFrame, N::Int)
    #get the vector of appropriate indices
    thermal_names = [string(thermal_builds[!,"Dim1"][i])*"_"*string(thermal_builds[!,"Dim2"][i]) for i=1:DataFrames.nrow(thermal_builds)]
    thermal_capacities = thermal_builds[!,"Val"]#want capacities

    #get regions and convert them to appropriate data where relevant
    thermal_categories = [string(thermal_builds[!,"Dim1"][i]) for i=1:DataFrames.nrow(thermal_builds)]
    thermal_regions = thermal_builds[!,"Dim2"]

    thermal_cap_factors = ones(length(thermal_names),N)

    return(create_generators_from_data(thermal_capacities.*thermal_cap_factors, thermal_names, thermal_categories, thermal_regions))
end

function process_vg(vg_builds::DataFrames.DataFrame,ReEDSfilepath::String,N::Int)
    #get the vector of appropriate indices
    vg_names = [string(vg_builds[!,"Dim1"][i])*"_"*string(vg_builds[!,"Dim2"][i]) for i=1:DataFrames.nrow(vg_builds)]
    vg_capacities = vg_builds[!,"Val"]#want capacities

    #get regions and convert them to appropriate data where relevant
    vg_categories = [string(vg_builds[!,"Dim1"][i]) for i=1:DataFrames.nrow(vg_builds)]
    vg_regions = vg_builds[!,"Dim2"]
    #pull in mapping info

    #once mapped, do we need to groupby/sum the resulting profiles as mapped?

    #load in recf
    cf_info = HDF5.h5read(joinpath(ReEDSfilepath,"inputs_case","recf.h5"), "data")
    indices = findfirst.(isequal.(vg_names), (cf_info["axis0"],))
    cap_factors = cf_info["block0_values"]
    # indices = findall(i->(i==(slicer-1)),load_info["axis1_label0"]) #I think b/c julia indexes from 1 we need -1 here
    retained_cap_factors = cap_factors[indices,1:N] #slice recf, just take first year for now until more is known
    
    #return the profiles
    return(create_generators_from_data(vg_capacities.*retained_cap_factors, vg_names, vg_categories, vg_regions))
end

function make_pras_system_from_mapping_info(ReEDSfilepath::String, Year::Int64)
    #######################################################
    # Loading the necessary mapping files and data
    #######################################################
    # Load the mapping metadata JSON file
    @info "Fetching ReEDS case data to build PRAS System..."

    #load-related information
    load_info = HDF5.h5read(joinpath(ReEDSfilepath,"inputs_case","load.h5"), "data") ;
    load_data = load_info["block0_values"];
    regions = load_info["axis0"];
    
    slicer = findfirst(isequal(Year), load_info["axis1_level0"]);
    indices = findall(i->(i==(slicer-1)),load_info["axis1_label0"]); #I think b/c julia indexes from 1 we need -1 here
    #probably want some kind of assert/error here if indices are empty that states the year is invalid
    load_year = load_data[:,indices]; #should be regionsX8760 if done right
    
    @info "Processing Areas in the mapping file into PRAS regions..."
    # area_names = string.(mapping_data[!,"PlainInter"]);
    num_areas = length(regions); 
    N = size(load_year)[2];#DataFrames.DataFrames.nrow(load_data)
    region_load = Array{Int64,2}(undef,num_areas,N);

    for (idx,region) in enumerate(regions)
        region_load[idx,:]=floor.(Int,load_year[idx,:]); #converts to int
    end
    new_regions = PRAS.Regions{N,PRAS.MW}(regions, region_load);
    
    #lines
    line_base_cap_data = process_lines(ReEDSfilepath);

    #generation capacity
    @info "splitting thermal, storage, vg generator types from installed ReEDS capacities..."
    thermal,storage,vg = split_generator_types(ReEDSfilepath,Year);
    
    #pull the generators into appropriate STRUCTS, with defaults

    #for thermal, we need outage stuff
    #for thermal, we will need a disaggreggation helper fxn. Possibly a couple. May need EIA860 data
    #for vg, we need profiles
    vg_tup = process_vg(vg,ReEDSfilepath,N);
    thermal_tup = process_thermals(thermal,N);

    gen_tup = [vcat(a,b) for (a,b) in zip(vg_tup,thermal_tup)];
    new_generators = PRAS.Generators{N,1,PRAS.Hour,PRAS.MW}(gen_tup[1],gen_tup[2],gen_tup[3],gen_tup[4],gen_tup[5]); #there's probably a way to better unpack this

    #######################################################
    # PRAS Timestamps
    #######################################################
    @info "Processing PRAS timestamps..."
    start_datetime = Dates.DateTime(string(Year)*"-01-01");
    finish_datetime = start_datetime + Dates.Hour(N-1);
    my_timestamps = StepRange(start_datetime, Dates.Hour(1), finish_datetime);

    #######################################################
    # PRAS Storages
    # No Storages in this system for now
    #######################################################
    @info "Processing Storages [EMPTY FOR NOW]..."
    stor_names = String[];
    stor_categories = String[];

    n_stor = 0;

    stor_charge_cap_array = Matrix{Int64}(undef, n_stor, N);
    stor_discharge_cap_array = Matrix{Int64}(undef, n_stor, N);
    stor_energy_cap_array = Matrix{Int64}(undef, n_stor, N);
    stor_chrg_eff_array = Matrix{Float64}(undef, n_stor, N);
    stor_dischrg_eff_array  = Matrix{Float64}(undef, n_stor, N);
    stor_cryovr_eff = ones(n_stor,N);
    λ_stor = Matrix{Float64}(undef, n_stor, N); 
    μ_stor = Matrix{Float64}(undef, n_stor, N);

    new_storage = PRAS.Storages{N,1,PRAS.Hour,PRAS.MW,PRAS.MWh}(stor_names,stor_categories,
                                                stor_charge_cap_array,stor_discharge_cap_array,stor_energy_cap_array,
                                                stor_chrg_eff_array,stor_dischrg_eff_array, stor_cryovr_eff,
                                                λ_stor,μ_stor);

    area_stor_idxs = fill(1:0, num_areas);

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

    # pras_system = PRAS.SystemModel(new_regions, new_interfaces, new_generators, area_gen_idxs, new_storage, area_stor_idxs, new_gen_stors,
    #                         area_genstor_idxs, new_lines,interface_line_idxs,my_timestamps);

    # return(pras_system)
    return(new_gen_stors)
end