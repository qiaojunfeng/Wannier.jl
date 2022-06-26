#!/usr/bin/env julia
using Wannier

function main()

    # seed_name = "silicon"
    seed_name = "/home/jqiao/git/Wannier.jl/test/fixtures/silicon"

    model = read_seedname(seed_name)

    A = max_localize(model)

    write_amn("silicon.maxloc.amn", A)
end


main()
