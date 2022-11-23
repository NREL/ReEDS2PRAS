
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

    vg_types = clean_names(vg_types)
    storage_types = clean_names(storage_types)
    vg_types = expand_types(vg_types,15)

    capacities = joinpath(ReEDSfilepath,"outputs","cap.csv")
    capacity_data = DataFrames.DataFrame(CSV.File(capacities))
    annual_data = capacity_data[capacity_data.Dim3.==Year,:] #only built capacity in the PRAS-year, since ReEDS capacity is cumulative

    #clear vg capacity on a regex, though there might be a better way...

    vg_capacity = annual_data[(findall(in(vg_types),annual_data.Dim1)),:]
    storage_capacity = annual_data[(findall(in(storage_types),annual_data.Dim1)),:]
    thermal_capacity = annual_data[(findall(!in(vcat(vg_types,storage_types)),annual_data.Dim1)),:]

    return(thermal_capacity,storage_capacity,vg_capacity)
end

function make_pras_system_from_mapping_info(ReEDSfilepath::String, Year::Int64)
    #######################################################
    # Loading the necessary mapping files and data
    #######################################################
    # Load the mapping metadata JSON file
    @info "Fetching ReEDS case data to build PRAS System..."

    #load-related information
    load_info = HDF5.h5read(joinpath(ReEDSfilepath,"inputs_case","load.h5"), "data") 
    load_data = load_info["block0_values"]
    regions = load_info["axis0"]
    
    slicer = findfirst(isequal(Year), load_info["axis1_level0"])
    indices = findall(i->(i==(slicer-1)),load_info["axis1_label0"]) #I think b/c julia indexes from 1 we need -1 here
    #probably want some kind of assert/error here if indices are empty that states the year is invalid
    load_year = load_data[:,indices] #should be regionsX8760 if done right
    # println(typeof(load_year))
    # new_regions = PRAS.Regions{N,PRAS.MW}(regions, load_year);
    
    @info "Processing Areas in the mapping file..."
    # area_names = string.(mapping_data[!,"PlainInter"]);
    num_areas = length(regions); 
    N = size(load_year)[2]#DataFrames.nrow(load_data)
    region_load = Array{Int64,2}(undef,num_areas,N);

    for (idx,region) in enumerate(regions)
        region_load[idx,:]=floor.(Int,load_year[idx,:]); #converts to int
    end
    new_regions = PRAS.Regions{N,PRAS.MW}(regions, region_load);
    
    #lines
    line_base_cap_data = process_lines(ReEDSfilepath)

    #generation capacity
    @info "splitting thermal, storage, vg generator types..."
    split_generator_types(ReEDSfilepath,Year)

    #we will need a disaggreggation helper fxn. Possibly a couple
    #possibly need to load in EIA860 data?

    #for thermal, we need outage stuff

    #for vg, we need profiles
    return(line_base_cap_data)
    # return (regions,slicer,indices,size(load_year))
end