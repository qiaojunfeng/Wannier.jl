"""
Generate parallel transport gauge Wannier functions.

# Args

- `seedname`: seedname for WIN/AMN/MMN/EIG files

# Options

- `-o, --output`: filename for output AMN. Default is `seedname.ptg.amn`
"""
@cast function ptg(seedname::String; output::String="")
    if output == ""
        output = basename(seedname) * ".ptg.amn"
    end

    model = read_w90(seedname; amn=false, eig=false)

    n_bands = model.n_bands
    n_bvecs = model.n_bvecs
    n_k1, n_k2, n_k3 = model.kgrid
    @info "Parallel transport gauge" n_bands n_k1 n_k2 n_k3 n_bvecs

    #Build the gauge that makes the Bloch frame continuous on the Brillouin Zone.
    #This is equivalent to building a set of algebraic decaying Wannier functions
    A, obs = parallel_transport(model; log_interp=false)

    Ωⁱ = omega(model.bvectors, model.M, model.A)
    @info "Initial spread"
    show(Ωⁱ)
    println("\n")

    Ωᶠ = omega(model.bvectors, model.M, A)
    @info "Final spread"
    show(Ωᶠ)
    println("\n")

    # Plot and display the results
    # Wannier.plot_obstruction(model, A, obs)

    write_amn(output, A)

    return nothing
end
