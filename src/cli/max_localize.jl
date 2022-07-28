"""
Maximally localize an isolated group of bands.

# Args

- `seedname`: seedname for WIN/AMN/MMN/EIG files

# Options

- `-o, --output`: filename for output AMN. Default is `seedname.maxloc.amn`
"""
@cast function maxloc(seedname::String; output::String="")
    if output == ""
        output = basename(seedname) * ".maxloc.amn"
    end

    model = read_seedname(seedname)

    A = max_localize(model)

    write_amn(output, A)

    return nothing
end

# if abspath(PROGRAM_FILE) == @__FILE__
#     main()
# end
