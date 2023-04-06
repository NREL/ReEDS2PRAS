# Converting FOR and MTTR to λ and μ
"""
    This function calculates the outage rate of a generator based on the
    forced outage rate and its Mean Time To Repair (MTTR).

    Parameters
    ----------
    for_gen: Float64
        forced outage rate (amount of time unit is out-of-service)
    mttr: Int64
        Mean time to repair (in hours)

    Returns
    -------
    (λ, μ): Tuple[Float64, Float64]
        λ is the probability of a unit being down, μ is probability 
        of recovery
"""
function outage_to_rate(for_gen::Float64, mttr::Int64)
    if (for_gen == 0.0)
        return (λ = 0.0, μ = 1.0)
    else
        if (for_gen > 1.0)
            for_gen = for_gen / 100
        end

        if (mttr != 0)
            μ = 1 / mttr
        else
            μ = 0.0
        end
        λ = (μ * for_gen) / (1 - for_gen)

        return (λ = λ, μ = μ)
    end
end

emptyvec(::Vector{<:Storage}) = Storage[]

get_components(comps::Vector{<:Storage}, region_name::String) =
    get_storages_in_region(comps, region_name)

emptyvec(::Vector{<:Generator}) = Generator[]

get_components(comps::Vector{<:Generator}, region_name::String) =
    get_generators_in_region(comps, region_name)
# Functions for processing ReEDS2PRAS generators and storages to prepare
# them for PRAS System generation
"""
    Gets components in each region of the system and reorganizes them into
    single sorted component vector and corresponding component index vector.

    Parameters
    ----------
    comps : COMPONENTS
        Vector containing components for each region
    region_names : Vector{String}
        Vector with names of regions in the system

    Returns
    -------
    sorted_comps : COMPONENTS
        Sorted vector of all components from every region
    region_comp_idxs : UnitRange{Int64}, 1
        Index vector pointing to components belonging to each specified region
"""
function get_sorted_components(
    comps::COMPONENTS,
    region_names::Vector{String},
) where {COMPONENTS <: Union{Vector{<:Generator}, Vector{<:Storage}}}
    num_regions = length(region_names)
    all_comps = []
    start_idx = Array{Int64}(undef, num_regions)
    region_comp_idxs = Array{UnitRange{Int64}, 1}(undef, num_regions)

    for (idx, region_name) in enumerate(region_names)
        region_comps = get_components(comps, region_name)
        push!(all_comps, region_comps)
        if idx == 1
            start_idx[idx] = 1
        else
            prev_idx = start_idx[idx - 1]
            prev_length = length(all_comps[idx - 1])
            start_idx[idx] = prev_idx + prev_length
        end
        region_comp_idxs[idx] = range(start_idx[idx], length = length(all_comps[idx]))
    end

    sorted_comps = emptyvec(comps)
    for idx in eachindex(all_comps)
        if (length(all_comps[idx]) != 0)
            append!(sorted_comps, all_comps[idx])
        end
    end
    return sorted_comps, region_comp_idxs
end

get_sorted_components(
    comps::COMPONENTS,
    regions::Vector{Region},
) where {COMPONENTS <: Union{Vector{<:Generator}, Vector{<:Storage}}} =
    get_sorted_components(comps, get_name.(regions))


# Functions for  processing ReEDS2PRAS lines (preparing PRAS lines)

"""
    Returns a list of tuples sorted by region name.

    Parameters
    ----------
    lines : Vector{Line}
        A list of lines containing the regions from and to.

    Returns
    -------
    List[Tuple[str]]
        A list of tuples sorted by region name.
"""
function get_sorted_region_tuples(lines::Vector{Line}, region_names::Vector{String})
    region_idxs = Dict(name => idx for (idx, name) in enumerate(region_names))

    line_from_to_reg_idxs = similar(lines, Tuple{Int, Int})

    for (l, line) in enumerate(lines)
        from_name = get_region_from(line)
        to_name = get_region_to(line)

        from_idx = region_idxs[from_name]
        to_idx = region_idxs[to_name]

        line_from_to_reg_idxs[l] =
            from_idx < to_idx ? (from_idx, to_idx) : (to_idx, from_idx)
    end

    return line_from_to_reg_idxs
end

function get_sorted_region_tuples(lines::Vector{Line}, regions::Vector{Region})
    get_sorted_region_tuples(lines, get_name.(regions))
end

function get_sorted_region_tuples(lines::Vector{Line})
    regions_from = get_region_from.(lines)
    regions_to = get_region_to.(lines)

    region_names = unique(append!(regions_from, regions_to))

    get_sorted_region_tuples(lines, region_names)
end

"""
    Returns a list of lines sorted 
    (acc. to sorted regions tuples and interface_region_idxs and interface_line_idxs for PRAS)

    Parameters
    ----------
    lines : Vector{Line}
        A Vector of ReEDS2PRAS Line objects
    region_names : Vector{String}
        Vector with names of regions in the system

    Returns
    -------
    Vector{Line}, UnitRange{Int64,1}, UnitRange{Int64,1} 
        A list of sorted lines, Index vector pointing to interface region_from and to belonging to each specified region,
        Index vector pointing to Lines belonging to an interface
"""

function get_sorted_lines(lines::Vector{Line}, region_names::Vector{String})
    line_from_to_reg_idxs = get_sorted_region_tuples(lines, region_names)
    line_ordering = sortperm(line_from_to_reg_idxs)

    sorted_lines = lines[line_ordering]
    sorted_from_to_reg_idxs = line_from_to_reg_idxs[line_ordering]
    interface_reg_idxs = unique(sorted_from_to_reg_idxs)

    # Ref tells Julia to use interfaces as Vector, only broadcasting over
    # lines_sorted
    interface_line_idxs = searchsorted.(Ref(sorted_from_to_reg_idxs), interface_reg_idxs)

    return sorted_lines, interface_reg_idxs, interface_line_idxs
end

get_sorted_lines(lines::Vector{Line}, regions::Vector{Region}) =
    get_sorted_lines(lines, get_name.(regions))

    """
    This code takes in a vector of Lines and a vector of Regions as input
    parameters. It filters the Lines to find VSC (voltage source converter)
    lines and non-VSC lines. Then, for each VSC line, it creates two new Line
    objects representing direct current (DC) converter capacityfor the regional connections 
    to the VSC line. Finally, a vector of all lines and a vector of all 
    regions are returned as output.

    Parameters
    ----------
    lines : Vector[Line]
        Vector of Line objects that contain information about all the lines in
        the system
    regions : Vector[Region]
        Vector of Region objects that contain information about all regions in
        the system

    Returns
    -------
    non_vsc_dc_lines : Vector[Line]
        Vector of Line objects with non_vsc_dc_lines
    regions : Vector[Region]
        Vector of Region objects with added DC lines
"""
function process_vsc_lines(lines::Vector{Line}, regions::Vector{Region})
    timesteps = first(regions).timesteps
    non_vsc_dc_lines = filter(line -> ~line.VSC, lines)
    vsc_dc_lines = filter(line -> line.VSC, lines)

    for vsc_line in vsc_dc_lines
        dc_region_from = "DC_$(vsc_line.region_from)"
        dc_region_to = "DC_$(vsc_line.region_to)"

        for reg_name in [dc_region_from, dc_region_to]
            if ~(reg_name in get_name.(regions))
                push!(regions, Region(reg_name, timesteps, zeros(Float64, timesteps)))
            end
        end

        push!(
            non_vsc_dc_lines,
            Line(
                "$(vsc_line.name)_DC",
                vsc_line.timesteps,
                vsc_line.category,
                dc_region_from,
                dc_region_to,
                vsc_line.forward_cap,
                vsc_line.backward_cap,
                vsc_line.legacy,
                vsc_line.FOR,
                vsc_line.MTTR,
            ),
        )
        push!(
            non_vsc_dc_lines,
            Line(
                "$(vsc_line.name)_Converter_From",
                vsc_line.timesteps,
                vsc_line.category,
                dc_region_from,
                vsc_line.region_from,
                vsc_line.converter_capacity[vsc_line.region_from],
                vsc_line.converter_capacity[vsc_line.region_from],
                vsc_line.legacy,
                vsc_line.FOR,
                vsc_line.MTTR,
            ),
        )
        push!(
            non_vsc_dc_lines,
            Line(
                "$(vsc_line.name)_Converter_To",
                vsc_line.timesteps,
                vsc_line.category,
                dc_region_to,
                vsc_line.region_to,
                vsc_line.converter_capacity[vsc_line.region_to],
                vsc_line.converter_capacity[vsc_line.region_to],
                vsc_line.legacy,
                vsc_line.FOR,
                vsc_line.MTTR,
            ),
        )
    end
    return non_vsc_dc_lines, regions
end
