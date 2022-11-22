
#####################################################
# Luke
# NREL
# November 2022
# ReEDS2PRAS for NTP
# Make a PRAS System from ReEDS H5s and CSVs

function make_pras_system_from_mapping_info(ReEDSfilepath::String, Year::Int64)
    #######################################################
    # Loading the necessary mapping files and data
    #######################################################
    # Load the mapping metadata JSON file
    @info "Fetching ReEDS case data to build PRAS System..."

    #load-related information
    load_info = h5read(joinpath(ReEDSfilepath,"load.h5"), "data") 
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
    line_cap_base_csv_location = joinpath(ReEDSfilepath,"trancap_init_energy.csv") #there is also a _prm file, not sure which is right?
    line_cap_additional_csv_location = joinpath(ReEDSfilepath,"trancap_fut.csv")

    line_base_cap_data = DataFrames.DataFrame(CSV.File(line_cap_base_csv_location;select=["r","rr","trtype","MW"]));
    line_additional_cap_data = DataFrames.DataFrame(CSV.File(line_cap_additional_csv_location;select=["r","rr","trtype","MW"]));
    
    #generation capacity

    #we will need a disaggreggation helper fxn. Possibly a couple
    #possibly need to load in EIA860 data?

    #for thermal, we need outage stuff

    #for vg, we need profiles
    return(line_base_cap_data,line_additional_cap_data)
    # return (regions,slicer,indices,size(load_year))
end