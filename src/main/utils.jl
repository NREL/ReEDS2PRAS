function clean_names!(input_vec::Vector{<:AbstractString})
    for (idx,a) in enumerate(input_vec)
        if occursin("*",a)
            input_vec[idx] = match(r"\*([a-zA-Z]+-*[a-zA-Z]*)_*", a)[1]
        end
    end 
    return input_vec
end

function expand_types!(input_vec::Vector{<:AbstractString},N::Int64)
    add_names = vec(["$(a)_$(n)" for a in input_vec, n in 1:N])
    return vcat(input_vec,add_names)
end

function run_pras_system(sys::PRAS.SystemModel,sample::Int)
    shortfalls,flows = PRAS.assess(sys,PRAS.SequentialMonteCarlo(samples=sample,seed=1),PRAS.Shortfall(),PRAS.Flow())
    println(PRAS.LOLE(shortfalls))
    println(PRAS.EUE(shortfalls))
    return shortfalls,flows
end

function Load_EIA_NEMS_DB(ReEDS_directory::String)
    EIA_NEMS_loc = joinpath(ReEDS_directory,"inputs","capacitydata","ReEDS_generator_database_final_EIA-NEMS.csv"); #there is also a _prm file, not sure which is right?
    EIA_NEMS_data = DataFrames.DataFrame(CSV.File(EIA_NEMS_loc));
    return EIA_NEMS_data
end

function disagg_existing_capacity(eia_df::DataFrames.DataFrame,built_capacity::Int,tech::String,pca::String,gen_for::Float64,N::Int,Year::Int)
    
    MTTR = 24;
    # @info "capacity in for $tech $pca is $built_capacity MW"
    tech_ba_year_existing = eia_df[(eia_df.tech.==tech) .& (eia_df.reeds_ba.==pca) .& (eia_df.RetireYear.>=Year) .& (eia_df.StartYear.<=Year), :];
    
    if length(tech_ba_year_existing.tech)==0
        return [thermal_gen("$(tech)_$(pca)_1",N,pca,built_capacity,tech,"New",gen_for,MTTR)]
    end

    remaining_capacity = built_capacity;
    existing_capacity = tech_ba_year_existing[!,"cap"]
    
    tech_len = length(existing_capacity);
    max_cap = maximum(existing_capacity)
    avg_cap = Statistics.mean(existing_capacity)

    generators_array = [];
    for (idx,built_cap) in enumerate(existing_capacity)
        int_built_cap = floor.(Int,built_cap);
        if int_built_cap < remaining_capacity
            remaining_capacity = remaining_capacity - int_built_cap;
            gen = thermal_gen("$(tech)_$(pca)_$(idx)",N,pca,int_built_cap,tech,"Existing",gen_for,MTTR)
            push!(generators_array,gen);
        else
            gen = thermal_gen("$(tech)_$(pca)_$(idx)",N,pca,remaining_capacity,tech,"Existing",gen_for,MTTR);
            push!(generators_array,gen);
            remaining_capacity = 0;
            break
        end 
    end
    #whatever remains, we want to build as new capacity

    if remaining_capacity > 0
        generators_array = disagg_new_capacity(generators_array,remaining_capacity,floor.(Int,avg_cap),floor.(Int,max_cap),tech,pca,gen_for,N,Year,MTTR);
    end

    return generators_array
end

function disagg_new_capacity(generators_array::Vector{<:Any},new_capacity::Int,avg::Int,max::Int,tech::String,pca::String,gen_for::Float64,N::Int,Year::Int,MTTR::Int)
    if avg==0
        return push!(generators_array,thermal_gen("$(tech)_$(pca)_new_1",N,pca,new_capacity,tech,"New",gen_for,MTTR))
    end
    n_gens = floor.(Int,new_capacity/avg);
    if n_gens==0
        return push!(generators_array,thermal_gen("$(tech)_$(pca)_new_1",N,pca,new_capacity,tech,"New",gen_for,MTTR))
    end
    remainder = new_capacity-(n_gens*avg);
    addtl_cap_per_gen = floor.(Int,remainder/n_gens);
    per_gen_cap = avg+addtl_cap_per_gen;
    small_remainder = new_capacity-(n_gens*per_gen_cap)

    for i in range(1,n_gens)
        push!(generators_array,thermal_gen("$(tech)_$(pca)_new_$(i)",N,pca,per_gen_cap,tech,"New",gen_for,MTTR))
    end

    push!(generators_array,thermal_gen("$(tech)_$(pca)_new_$(n_gens+1)",N,pca,small_remainder,tech,"New",gen_for,MTTR)); #integer remainder is made into a tiny gen
    return generators_array
end

abstract type CEMdata end
struct ReEDSdata <:CEMdata
    ReEDSfilepath::String
    Year::Int

    # Inner Constructors
    # ReEDSdata(nothing) = ReEDSdata("/projects/ntps/llavin/ReEDS-2.0/runs/erc_conv_2_ercot_seq",2028)
    # arguably I should get rid of this - 
    # Checks
    ReEDSdata(x,y) =
    if ~(2020 < y <= 2050)
        error("Year should be between 2020 and 2050 for ReEDS case for now")
    else
        new(x,y)
    end
end
function get_load_file(data::ReEDSdata)
    if !isfile(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","plot_load_"*string(data.Year)*".h5"))
        Year = data.Year;
        error("The year $Year does not have an associated Augur load h5 file. Are you sure ReeDS was run and Augur results saved for $Year?")
    end
    return HDF5.h5read(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","plot_load_"*string(data.Year)*".h5"),"data")
end

function get_vg_cf_data(data::ReEDSdata)
    if !isfile(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","plot_vre_gen_"*string(data.Year)*".h5"))
        Year = data.Year;
        error("The year $Year does not have an associated Augur vg h5 file. Are you sure ReeDS was run and Augur results saved for $Year?")
    end
    return HDF5.h5read(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","plot_vre_gen_"*string(data.Year)*".h5"),"data")
end

function get_forced_outage_data(data::ReEDSdata)
    if !isfile(joinpath(data.ReEDSfilepath,"inputs_case","outage_forced.csv"))
        error("No forced outage data is found.")
    end    
    return DataFrames.DataFrame(CSV.File(joinpath(data.ReEDSfilepath,"inputs_case","outage_forced.csv"),header=false))
end

function get_line_capacity_data(data::ReEDSdata)
    if !isfile(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","tran_cap_"*string(data.Year)*".csv"))
        Year = data.Year;
        error("The year $Year does not have transmission capacity data. Are you sure ReEDS was run and Augur results saved for $Year?")
    end
    return DataFrames.DataFrame(CSV.File(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","tran_cap_"*string(data.Year)*".csv")));#load the csv
end

function get_converter_capacity_data(data::ReEDSdata)
    if !isfile(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","cap_converter_"*string(data.Year)*".csv"))
        Year = data.Year;
        error("The year $Year does not have capacity converter data. Are you sure ReEDS was run and Augur results saved for $Year?")
    end
    return DataFrames.DataFrame(CSV.File(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","cap_converter_"*string(data.Year)*".csv")));
end


function get_technology_types(data::ReEDSdata)
    if !isfile(joinpath(data.ReEDSfilepath,"inputs_case","tech-subset-table.csv"))
        error("no table of technology types!")
    end
    return DataFrames.DataFrame(CSV.File(joinpath(data.ReEDSfilepath,"inputs_case","tech-subset-table.csv")))
end

function get_region_mapping(data::ReEDSdata)
    if !isfile(joinpath(data.ReEDSfilepath,"inputs_case","rsmap.csv"))
        error("no table of r-s region mapping!")
    end
    return DataFrames.DataFrame(CSV.File(joinpath(data.ReEDSfilepath,"inputs_case","rsmap.csv")))
end

function get_ICAP_data(data::ReEDSdata)
    if !isfile(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","max_cap_"*string(data.Year)*".csv"))
        Year = data.Year;
        error("The year $Year does not have generator installed capacity data. Are you sure REEDS was run and Augur results saved for year $Year")
    end
    return DataFrames.DataFrame(CSV.File(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","max_cap_"*string(data.Year)*".csv")))
end

function get_storage_energy_capacity_data(data::ReEDSdata)
    if !isfile(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","energy_cap_"*string(data.Year)*".csv"))
        Year = data.Year;
        error("The year $Year does not have generator installed storage energy capacity data. Are you sure REEDS was run and Augur results saved for year $Year")
    end
    return DataFrames.DataFrame(CSV.File(joinpath(data.ReEDSfilepath,"ReEDS_Augur","augur_data","energy_cap_"*string(data.Year)*".csv")))
end