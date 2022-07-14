"""
Generate parallel transport gauge Wannier functions.

# Args

- `seedname`: seedname for WIN/AMN/MMN/EIG files

# Options

- `-o, --output`: filename for output AMN. Default is `seedname.ptg.amn`
"""
@cast function ptg(seedname::String; output::Union{String,Nothing}=nothing)

    # seedname = "tests/KaneMele/kanemele_0.0_25"
    # seedname = "tests/silicon/silicon"
    # seedname = "/home/jqiao/git/Wannier.jl/test/fixtures/val/silicon"
    if output === nothing
        output = basename(seedname) * ".ptg.amn"
    end

    model = read_seedname(seedname; amn=false, eig=false)

    n_bands = model.n_bands
    n_bvecs = model.n_bvecs
    n_k1, n_k2, n_k3 = model.kgrid
    @info "Parallel transport gauge" n_bands n_k1 n_k2 n_k3 n_bvecs

    #Build the gauge that makes the Bloch frame continuous on the Brillouin Zone.
    #This is equivalent to building a set of algebraic decaying Wannier functions
    A, obs = parallel_transport(model; log_interp=false)

    Ωⁱ = omega(model.bvectors, model.M, model.A)
    @info "Initial spread"
    pprint(Ωⁱ)

    Ωᶠ = omega(model.bvectors, model.M, A)
    @info "Final spread"
    pprint(Ωᶠ)

    # Plot and display the results
    # Wannier.plot_obstruction(model, A, obs)

    write_amn(output, A)

    return nothing
end
