"""
Maximally localize an isolated group of bands.

# Args

- `seedname`: seedname for `win`/`amn`/`mmn`/`eig` files

# Options

- `-o, --output=<str>`: filename for output `amn`. Default is `seedname.maxloc.amn`
- `-m, --maxiter=<int>`: max number of iterations. Default is `50`
"""
@cast function maxloc(seedname::String; output::String="", maxiter::Int=50)
    if output == ""
        output = basename(seedname) * ".maxloc.amn"
    end

    model = read_w90(seedname)

    A = max_localize(model; max_iter=maxiter)

    write_amn(output, A)

    return nothing
end

# if abspath(PROGRAM_FILE) == @__FILE__
#     main()
# end
