# Run following from the repo base directory
main_paths = ["./src", "./test"]
for main_path in main_paths
    for (root, dir, files) in walkdir(main_path)
        for f in files
        @show file_path = abspath(root, f)
        !occursin(".jl", f) && continue
        occursin("generated", file_path) && continue
        format(file_path;
            whitespace_ops_in_indices = true,
            remove_extra_newlines = true,
            verbose = true,
            always_for_in = true,
            whitespace_typedefs = true
            )
        end
    end
end