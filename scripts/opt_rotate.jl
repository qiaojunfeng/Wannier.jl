#!/usr/bin/env julia
using Wannier

function main()

    # seedname = "silicon"
    seedname = "/home/jqiao/git/Wannier.jl/test/test_grad/silicon"

    # Input AMN is parallel transport gauge
    model = read_seedname(seedname)

    A = opt_rotate(model)

    write_amn("silicon.optrot.amn", A)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
