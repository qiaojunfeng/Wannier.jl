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
    @info "Initial total spread" Ω = round(Ωⁱ.Ω; digits = 5)
    @info "Initial spread:" ω = round.(Ωⁱ.ω'; digits = 5)
    @info "Initial centers:" r = round.(Ωⁱ.r; digits = 5)

    Ωᶠ = omega(model.bvectors, model.M, A)
    @info "Final total spread" Ω = round(Ωᶠ.Ω; digits = 5)
    @info "Final spread:" ω = round.(Ωᶠ.ω'; digits = 5)
    @info "Final centers:" r = round.(Ωᶠ.r; digits = 5)

    # Plot and display the results
    Wannier.plot_obstruction(model, A, obs)

    write_amn("silicon.ptg.amn", A)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
