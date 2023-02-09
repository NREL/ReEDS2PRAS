function load_objects(ReEDS_data,WEATHERYEAR::Int,N::Int,year::Int,NEMS_path::String,min_year::Int)
    
    @info "Processing regions and associating load profiles..."
    region_array = process_regions_and_load(ReEDS_data,WEATHERYEAR,N)
    
    @info "Processing lines and adding VSC-related regions, if applicable..."
    lines = process_lines(ReEDS_data,get_name.(region_array),year,N)
    lines,regions = process_vsc_lines(lines,region_array)
    
    # Create Generator Objects
    # **TODO: Should 0 MW generators be allowed after disaggregation?
    # **TODO: Should hydro be split out as a generator-storage?
    # **TODO: is it important to also handle planned outages?
    @info "splitting thermal, storage, vg generator types from installed ReEDS capacities..."
    thermal,storage,variable_gens = split_generator_types(ReEDS_data,year)

    @info "reading in ReEDS generator-type forced outage data..."
    forced_outage_data = get_forced_outage_data(ReEDS_data)
    forced_outage_dict = Dict(forced_outage_data[!,"ResourceType"] .=> forced_outage_data[!,"FOR"])

    @info "Processing conventional/thermal generators..."
    thermal_gens = process_thermals_with_disaggregation(ReEDS_data,thermal,forced_outage_dict,N,year,NEMS_path)
    @info "Processing variable generation..."
    gens_array = process_vg(thermal_gens,variable_gens,forced_outage_dict,ReEDS_data,year,WEATHERYEAR,N,min_year)

    @info "Processing Storages..."
    storage_array = process_storages(storage,forced_outage_dict,ReEDS_data,N,get_name.(regions),year)

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
    storages, stor_idxs = get_sorted_components(storages,regions)

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
    #**TODO: can get a bug/error here if ReEDS case lacks multiple load years
    slicer = findfirst(isequal(weather_year), load_info["axis1_level1"]) #slices based on input weather year
    indices = findall(i->(i==(slicer-1)),load_info["axis1_label1"]) #I think b/c julia indexes from 1 we need -1 here
    
    load_year = load_data[:,indices[1:N]] #should be regionsX8760, which is now just enforced

    return [Region(r,N,floor.(Int,load_year[idx,:])) for (idx,r) in enumerate(regions)]
end

function process_lines(ReEDS_data,regions::Vector{<:AbstractString},year::Int,N::Int)

    line_base_cap_data = get_line_capacity_data(ReEDS_data) #it is assumed this has prm line capacity data

    converter_capacity_data = get_converter_capacity_data(ReEDS_data)
    converter_capacity_dict = Dict(converter_capacity_data[!,"r"] .=> converter_capacity_data[!,"MW"])

    #add 0 converter capacity for regions that lack a converter
    if length(keys(converter_capacity_dict))>0
        for reg in regions
            if !(reg in keys(converter_capacity_dict))
                @info "$reg does not have VSC converter capacity, so adding a 0"
                converter_capacity_dict[String(reg)] = 0
            end
        end
    end

    function keep_line(from_pca, to_pca)
        from_idx = findfirst(x->x==from_pca,regions)
        to_idx = findfirst(x->x==to_pca,regions)
        return ~isnothing(from_idx) && ~isnothing(to_idx) && (from_idx < to_idx)
    end
    system_line_naming_data = DataFrames.subset(line_base_cap_data, [:r, :rr] => DataFrames.ByRow(keep_line))
    system_line_naming_data = DataFrames.combine(DataFrames.groupby(system_line_naming_data, ["r","rr","trtype"]), :MW => sum) #split-apply-combine b/c some lines have same name convention

    lines_array = Line[]
    for row in eachrow(system_line_naming_data)
        forward_cap = sum(line_base_cap_data[(line_base_cap_data.r.==row.r) .& (line_base_cap_data.rr.==row.rr) .& (line_base_cap_data.trtype.==row.trtype),"MW"])
        backward_cap = sum(line_base_cap_data[(line_base_cap_data.r.==row.rr) .& (line_base_cap_data.rr.==row.r) .& (line_base_cap_data.trtype.==row.trtype),"MW"])

        name = "$(row.r)_$(row.rr)_$(row.trtype)"        
        @info "a line $name, with $forward_cap MW fwd and $backward_cap bckwrd in $(row.trtype)"
        if row.trtype!="VSC"
            push!(lines_array,Line(name,N,row.trtype,row.r,row.rr,forward_cap,backward_cap,"Existing",0.0,24))
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
    vg_types = expand_vg_types!(vg_types,unique(vg_resource_types.i)) 

    vg_capacity = capacity_data[(findall(in(vg_types),capacity_data.i)),:]
    storage_capacity = capacity_data[(findall(in(storage_types),capacity_data.i)),:]
    thermal_capacity = capacity_data[(findall(!in(vcat(vg_types,storage_types)),capacity_data.i)),:]
    return(thermal_capacity,storage_capacity,vg_capacity)
end

function process_thermals_with_disaggregation(ReEDS_data,thermal_builds::DataFrames.DataFrame,FOR_dict::Dict,N::Int,year::Int,NEMS_path::String)#FOR_data::DataFrames.DataFrame,
    thermal_builds = thermal_builds[(thermal_builds.i.!= "csp-ns"), :] #csp-ns is not a thermal; just drop in for now
    thermal_builds = DataFrames.combine(DataFrames.groupby(thermal_builds, ["i","r"]), :MW => sum) #split-apply-combine to handle differently vintaged entries
    EIA_db = Load_EIA_NEMS_DB(NEMS_path)

    all_generators = Generator[]
    for row in eachrow(thermal_builds) #this loop gets the FOR for each build/tech

        if row.i in keys(FOR_dict)
            gen_for = FOR_dict[row.i]
        else
            gen_for = 0.00 #assume as 0 for gens dropped from ReEDS table
            @info "CONVENTIONAL GENERATION: for $(row.i), and region $(row.r), no gen_for is found in ReEDS forced outage data, so $gen_for is used"
        end

        generator_array = disagg_existing_capacity(EIA_db,floor(Int,row.MW_sum),String(row.i),String(row.r),gen_for,N,year)
        append!(all_generators,generator_array)
    end
    return all_generators
end

function process_vg(generators_array::Vector{<:ReEDS2PRAS.Generator},vg_builds::DataFrames.DataFrame,FOR_dict::Dict,ReEDS_data,year::Int,weather_year::Int,N::Int,min_year::Int)

    region_mapper_df = get_region_mapping(ReEDS_data)
    region_mapper_dict = Dict(region_mapper_df[!,"rs"] .=> region_mapper_df[!,"*r"])
    cf_info = get_vg_cf_data(ReEDS_data)
    
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
            gen_for = .00 #make this 0 for vg if no match
        end
        name = "$(name)_"
        
        push!(generators_array,VG_Gen(name,N,region,maximum(profile),profile,category,"New",gen_for,24)) 
    end
    return generators_array
end

function process_storages(storage_builds::DataFrames.DataFrame,FOR_dict::Dict,ReEDS_data,N::Int,regions::Vector{<:AbstractString},year::Int64)
    
    storage_energy_capacity_data = get_storage_energy_capacity_data(ReEDS_data)
    energy_capacity_df = DataFrames.combine(DataFrames.groupby(storage_energy_capacity_data, ["i","r"]), :MWh => sum) #split-apply-combine to handle differently vintaged entries

    storages_array = Storage[]
    for (idx,row) in enumerate(eachrow(storage_builds))
        name = "$(string(row.i))_$(string(row.r))"
        if string(row.i) in keys(FOR_dict)
            gen_for = FOR_dict[string(row.i)]
        else
            gen_for = 0.0
            @info "STORAGE: did not find FOR for storage $name $(row.r) $(row.i), so setting FOR to default value $gen_for"
        end 
        name = "$(name)_"#append for later matching

        int_duration = round(energy_capacity_df[idx,"MWh_sum"]/row.MW) #as per discussion w/ patrick, find duration of storage, then make energy capacity on that duration?
        push!(storages_array,Battery(name,N,string(row.r),string(row.i),row.MW,row.MW,round(Int,row.MW)*int_duration,"New",1,1,1,gen_for,24))
    end
    return storages_array
end

function process_genstors(regions::Vector{<:AbstractString},N::Int)
    
    gen_stors = Gen_Storage[Gen_Storage("blank_1", N, regions[1], "blank_genstor", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)] #empty for now
    
    return gen_stors
end

function disagg_existing_capacity(eia_df::DataFrames.DataFrame,built_capacity::Int,tech::String,pca::String,gen_for::Float64,N::Int,year::Int)
    
    MTTR = 24
    tech_ba_year_existing = DataFrames.subset(eia_df, :tech=>DataFrames.ByRow(==(tech)), :reeds_ba=>DataFrames.ByRow(==(pca)),:RetireYear=>DataFrames.ByRow(>=(year)),:StartYear=>DataFrames.ByRow(<=(year)))

    if DataFrames.nrow(tech_ba_year_existing)==0
        return [Thermal_Gen("$(tech)_$(pca)_1",N,pca,built_capacity,tech,"New",gen_for,MTTR)]
    end

    remaining_capacity = built_capacity
    existing_capacity = tech_ba_year_existing[!,"cap"]
    
    tech_len = length(existing_capacity)
    max_cap = maximum(existing_capacity)
    avg_cap = Statistics.mean(existing_capacity)

    generators_array = []
    for (idx,built_cap) in enumerate(existing_capacity)
        int_built_cap = floor(Int,built_cap)
        if int_built_cap < remaining_capacity
            gen_cap = int_built_cap
            remaining_capacity -= int_built_cap
        else
            gen_cap = remaining_capacity
            remaining_capacity = 0
        end
        gen = Thermal_Gen("$(tech)_$(pca)_$(idx)",N,pca,gen_cap,tech,"Existing",gen_for,MTTR)
        push!(generators_array,gen)
    end

    #whatever remains, we want to build as new capacity
    if remaining_capacity > 0
        add_new_capacity!(generators_array,remaining_capacity,floor.(Int,avg_cap),floor.(Int,max_cap),tech,pca,gen_for,N,year,MTTR)
    end

    return generators_array
end

function add_new_capacity!(generators_array::Vector{<:Any},new_capacity::Int,existing_avg_unit_cap::Int,max::Int,tech::String,pca::String,gen_for::Float64,N::Int,year::Int,MTTR::Int)

    if existing_avg_unit_cap==0 #if there are no existing units to determine size of new unit(s), build all new capacity as a single generator
        return push!(generators_array,Thermal_Gen("$(tech)_$(pca)_new_1",N,pca,new_capacity,tech,"New",gen_for,MTTR))
    end
    
    n_gens = floor(Int,new_capacity/existing_avg_unit_cap)
    if n_gens==0
        return push!(generators_array,Thermal_Gen("$(tech)_$(pca)_new_1",N,pca,new_capacity,tech,"New",gen_for,MTTR))
    end

    remainder = new_capacity-(n_gens*existing_avg_unit_cap)
    addtl_cap_per_gen = floor(Int,remainder/n_gens)
    per_gen_cap = existing_avg_unit_cap+addtl_cap_per_gen
    for i in range(1,n_gens)
        push!(generators_array,Thermal_Gen("$(tech)_$(pca)_new_$(i)",N,pca,per_gen_cap,tech,"New",gen_for,MTTR))
    end

    small_remainder = new_capacity-(n_gens*per_gen_cap)
    if small_remainder>0
        push!(generators_array,Thermal_Gen("$(tech)_$(pca)_new_$(n_gens+1)",N,pca,small_remainder,tech,"New",gen_for,MTTR)) #integer remainder is made into a tiny gen
    end

    return generators_array
end

function make_pras_interfaces(sorted_lines::Vector{Line},interface_reg_idxs::Vector{Tuple{Int64, Int64}},interface_line_idxs::Vector{UnitRange{Int64}},
    regions::Vector{Region})
    make_pras_interfaces(sorted_lines,interface_reg_idxs,interface_line_idxs, get_name.(regions))
end

function make_pras_interfaces(sorted_lines::Vector{Line},interface_reg_idxs::Vector{Tuple{Int64, Int64}},interface_line_idxs::Vector{UnitRange{Int64}},
    region_names::Vector{String})

    num_interfaces = length(interface_reg_idxs);
    interface_regions_from = first.(interface_reg_idxs);
    interface_regions_to = last.(interface_reg_idxs);

    N = first(sorted_lines).N

    # Lines
    line_names = get_name.(sorted_lines)
    line_cats = get_category.(sorted_lines)
    line_forward_cap = reduce(vcat,get_forward_capacity.(sorted_lines))
    line_backward_cap = reduce(vcat,get_backward_capacity.(sorted_lines))
    line_λ = reduce(vcat,get_λ.(sorted_lines))
    line_μ = reduce(vcat,get_μ.(sorted_lines))

    new_lines = PRAS.Lines{N,1,PRAS.Hour,PRAS.MW}(line_names, line_cats, line_forward_cap, line_backward_cap, line_λ ,line_μ);

    interface_forward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);
    interface_backward_capacity_array = Matrix{Int64}(undef, num_interfaces, N);

    for i in 1:num_interfaces
        interface_forward_capacity_array[i,:] =  sum(line_forward_cap[interface_line_idxs[i],:],dims=1)
        interface_backward_capacity_array[i,:] =  sum(line_backward_cap[interface_line_idxs[i],:],dims=1)
    end

    new_interfaces = PRAS.Interfaces{N,PRAS.MW}(interface_regions_from, interface_regions_to, interface_forward_capacity_array, interface_backward_capacity_array);

    return new_lines, new_interfaces

end