#!/usr/bin/env julia
using Wannier

function main()

    # seedname = "silicon"
    seedname = "/home/jqiao/git/Wannier.jl/test/fixtures/silicon"

    model = read_seedname(seedname)

    A = max_localize(model)

    write_amn("silicon.maxloc.amn", A)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
