"""
Maximally localize w.r.t a single unitary rotation of all the kpoints.

Usually should start from parallel transport gauge `amn`, where the gauge
are already smoothened w.r.t. kpoints. However, there is sitll a global
unitary transformation freedom, which will be minimized by this
optimal rotation function.

# Args

- `seedname`: seedname for `win`/`amn`/`mmn`/`eig` files

# Options

- `-o, --output=<str>`: filename for output `amn`. Default is `seedname.optrot.amn`
- `-m, --maxiter=<int>`: max number of iterations. Default is `50`
"""
@cast function optrot(seedname::String; output::String="", maxiter::Int=50)
    if output == ""
        output = basename(seedname) * ".optrot.amn"
    end

    # Input AMN is parallel transport gauge
    model = read_w90(seedname)
    U0 = deepcopy(model.U)

    W = opt_rotate(model; max_iter=maxiter)

    U = merge_gauge(U0, W)
    write_amn(output, U)

    return nothing
end
