using TOML: parsefile

"""
Split valence and conduction Wannier functions.

Usually start from a Wannierization of valence+conduction bands.
Then this command split WFs into two independent groups.

# Args

- `seedname`: seedname for `win`/`amn`/`mmn`/`eig` files

# Options

- `--nval`: number of valence WFs. Default is `n_wann ÷ 2`
- `--outdir-val`: dirname for output valence `amn`/`mmn`/`eig`. Default is `val`
- `--outdir-cond`: dirname for output conduction `amn`/`mmn`/`eig`. Default is `cond`
- `--config`: config file for `splitvc` command, e.g.
    ```toml
    [groups]
    indices = [ [ 1, 2,], [ 3, 4, 5, 6,], [ 7, 8,], ]
    outdirs = [ "val_1", "val_2", "cond_3",]
    ```

# Flags

- `--run-disentangle`: read `amn` and run disentangle first, otherwise read `chk` to
    get unitary matrices from `n_bands` to `n_wann`
- `--run-optrot`: max localize w.r.t. single unitary matrix after parallel transport.
    Should further reduce the spread and much closer to the true max localization.
- `--run-maxloc`: run a final max localize w.r.t. all kpoints.
    Should reach the true max localization.
- `--rotate-unk`: generate `unk` files for valence and conduction, for plotting WFs
- `--binary`: write `amn`/`mmn`/`eig`/`unk` in Fortran binary format
"""
@cast function splitvc(
    seedname::String;
    nval::Int=0,
    outdir_val::String="val",
    outdir_cond::String="cond",
    config::String="",
    run_disentangle::Bool=false,
    run_optrot::Bool=false,
    run_maxloc::Bool=false,
    rotate_unk::Bool=false,
    binary::Bool=false,
)
    # Input AMN is Silicon s,p projection
    model = read_w90(seedname; amn=false)

    if run_disentangle
        # You can also use disentangle to get a good gauge from initial projection
        model.U .= read_amn_ortho("$seedname.amn")
        model.U .= disentangle(model)
    else
        # Get max localized gauge from chk file, 1st try text format, then binary
        if isfile("$seedname.chk.fmt")
            chk = read_chk("$seedname.chk.fmt")
        else
            chk = read_chk("$seedname.chk")
        end
        # We replace the initial projection by the "good" max loc gauge
        model.U .= get_U(chk)
    end

    if config == ""
        (nval == 0) && (nval = model.n_wann ÷ 2)
        @info "number of valence WFs = $nval"
        indices = [1:nval, (nval + 1):(model.n_wann)]
        outdirs = [outdir_val, outdir_cond]
    else
        @info "reading config file: $config"
        groups = parsefile(config)["groups"]
        indices = groups["indices"]
        outdirs = groups["outdirs"]
    end
    println("Model will be split into $(length(indices)) groups")
    for (i, idxs) in enumerate(indices)
        println("  Group $i:")
        println("    indices: $(idxs)")
        println("    outdir : $(outdirs[i])")
    end

    @info "Original model initial spread"
    show(omega(model))
    println("\n")

    models_Us = split_wannierize(model, indices)

    for (i, (m, U)) in enumerate(models_Us)
        @info "Group $i after parallel transport:"
        show(omega(m))
        println("\n")

        if run_optrot
            @info "Run optimal rotation"
            println()
            W = opt_rotate(m)
            m.U .= merge_gauge(m.U, W)
        end

        if run_maxloc
            @info "Run max localization"
            println()
            m.U .= max_localize(m)
        end

        # Write files
        outdir = joinpath(dirname(seedname), outdirs[i])
        !isdir(outdir) && mkdir(outdir)
        seedname_i = joinpath(outdir, basename(seedname))
        write_w90(seedname_i, m; binary=binary)
    end

    # UNK files for plotting WFs
    if rotate_unk
        dir = dirname(seedname)
        dir == "" && (dir = ".")

        outdirs = [joinpath(dir, odir) for odir in outdirs]
        split_unk(dir, [mU[2] for mU in models_Us], outdirs; binary=binary)
    end

    return nothing
end
