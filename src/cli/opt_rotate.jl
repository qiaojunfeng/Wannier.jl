"""
Maximally localize w.r.t a single unitary rotation of all the kpoints.

Usually should start from parallel transport gauge AMN, where the gauge
are already smoothened w.r.t. kpoints. However, there is sitll a global
unitary transformation freedom, which will be minimized by this
optimal rotation function.

# Args

- `seedname`: seedname for WIN/AMN/MMN/EIG files

# Options

- `-o, --output`: filename for output AMN. Default is `seedname.optrot.amn`
"""
@cast function optrot(seedname::String; output::Union{String,Nothing}=nothing)

    # seedname = "/home/jqiao/git/Wannier.jl/test/test_grad/silicon"
    if output === nothing
        output = basename(seedname) * ".optrot.amn"
    end

    # Input AMN is parallel transport gauge
    model = read_seedname(seedname)
    A0 = deepcopy(model.A)

    W = opt_rotate(model)

    A = rotate_A(A0, W)
    write_amn(output, A)

    return nothing
end
