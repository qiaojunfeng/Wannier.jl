#!/usr/bin/env julia
using Comonicon


"""
Maximally localize an isolated group of bands.

# Args

- `seedname`: seedname for WIN/AMN/MMN/EIG files

# Options

- `-o, --output`: filename for output AMN. Default is `seedname.maxloc.amn`
"""
@cast function maxloc(seedname::String; output::Union{String,Nothing} = nothing)

    # seedname = "/home/jqiao/git/Wannier.jl/test/fixtures/silicon"
    if output === nothing
        output = basename(seedname) * ".maxloc.amn"
    end

    model = read_seedname(seedname)

    A = max_localize(model)

    write_amn(output, A)
end


# if abspath(PROGRAM_FILE) == @__FILE__
#     main()
# end
