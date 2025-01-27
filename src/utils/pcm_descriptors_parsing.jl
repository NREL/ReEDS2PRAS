isjson = endswith(".json");

inner_keys = ["avg_capacity_MW","mean_time_to_repair","forced_outage_rate"]

function parse_pcm_default_descriptors(file::String)
    if ~(isjson(file))
        @error "PCM default input descriptors file provided is not in JSON format. Exiting..."
    end

    # pcm_default_input_dict = JSON.parsefile(file)
    # pcm_default_inputs = Dict()
    # for key in keys(pcm_default_input_dict)
    #     push!(pcm_default_inputs, key => inner_key)
    #     for inner_key in inner_keys
    #         push!(
    #             pcm_default_inputs[key],
    #             key => inner_key => pcm_default_input_dict[key][inner_key],
    #         )
    #     end
    # end
    return JSON.parsefile(file)
end