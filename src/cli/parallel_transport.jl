#!/usr/bin/env julia
using Wannier


function main()

    # seedname = "tests/KaneMele/kanemele_0.0_25"
    # seedname = "tests/silicon/silicon"

    seedname = "/home/jqiao/git/Wannier.jl/test/fixtures/val/silicon"

    model = read_seedname(seedname; amn = false, eig = false)

    n_bands = model.n_bands
    n_bvecs = model.n_bvecs
    n_k1, n_k2, n_k3 = model.kgrid
    @info "Parallel transport gauge" n_bands n_k1 n_k2 n_k3 n_bvecs

    #Build the gauge that makes the Bloch frame continuous on the Brillouin Zone.
    #This is equivalent to building a set of algebraic decaying Wannier functions
    A, obs = parallel_transport(model; log_interp = false, return_obs = true)

    Ωⁱ = omega(model.bvectors, model.M, model.A)
    @info "Initial spread"
    print_spread(Ωⁱ)

    Ωᶠ = omega(model.bvectors, model.M, A)
    @info "Final spread"
    print_spread(Ωᶠ)

    # Plot and display the results
    # Wannier.plot_obstruction(model, A, obs)

    write_amn("silicon.ptg.amn", A)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
