function clean_names(input_vec::Vector)
    for (idx,a) in enumerate(input_vec)
        if occursin("*",a)
            input_vec[idx] = match(r"\*([a-zA-Z]+-*[a-zA-Z]*)_*", a)[1]
        end
    end 
    return input_vec
end

function expand_types(input_vec::Vector,N::Int64)
    add_names = [];
    for a in input_vec
        for n in 1:N
            mystr = string(a*"_"*string(n))
            push!(add_names,mystr)
        end
    end
    return vcat(input_vec,add_names)
end