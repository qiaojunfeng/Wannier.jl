"""
Split valence and conduction Wannier functions.

Usually start from a Wannierization of valence+conduction bands.
Then this command split WFs into two independent groups.

# Args

- `seedname`: seedname for WIN/AMN/MMN/EIG files

# Options

- `--nval`: number of valence WFs. Default is `n_wann ÷ 2`
- `--outdir-val`: dirname for output valence AMN/MMN/EIG. Default is `val`
- `--outdir-cond`: dirname for output conduction AMN/MMN/EIG. Default is `cond`

# Flags

- `--run-disentangle`: run disentangle first, otherwise read CHK to
    get unitary matrices from `n_bands` to `n_wann`
- `--run-optrot`: max localize w.r.t. single unitary matrix after parallel transport.
    Should further reduce the spread and much closer to the true max localization.
- `--run-maxloc`: run a final max localize w.r.t. all kpoints.
    Should reach the true max localization.
- `--rotate-unk`: generate UNK files for valence and conduction, for plotting WFs
"""
@cast function splitvc(
    seedname::String;
    nval::Union{Int,Nothing}=nothing,
    outdir_val::String="val",
    outdir_cond::String="cond",
    run_disentangle::Bool=false,
    run_optrot::Bool=false,
    run_maxloc::Bool=false,
    rotate_unk::Bool=false,
)
    # seedname = "silicon"
    # seedname = "/home/jqiao/git/Wannier.jl/test/fixtures/silicon"

    # Input AMN is Silicon s,p projection
    model = read_seedname(seedname)

    # calculate spread
    f(m::Model) = omega(m.bvectors, m.M, m.A).Ω

    if run_disentangle
        # You can also use disentangle to get a good gauge from initial projection
        model.A = disentangle(model)
    else
        # Get max localized gauge from chk file
        chk = read_chk("$seedname.chk.fmt")
        # We replace the initial projection by the "good" max loc gauge
        model.A = get_amn(chk)
    end

    @info "Valence + conduction initial spread"
    pprint(f(model))

    (nval === nothing) && (nval = model.n_wann ÷ 2)

    # UNK files for plotting WFs
    splitted = split_wannierize(model, nval, rotate_unk)
    if rotate_unk
        model_val, model_cond, Uv, Uc = splitted
    else
        model_val, model_cond = splitted
    end

    @info "Valence after parallel transport:"
    pprint(f(model_val))
    @info "Conduction after parallel transport:"
    pprint(f(model_cond))

    if run_optrot
        @info "Run optimal rotation"
        println()
        Wv = opt_rotate(model_val)
        model_val.A .= rotate_amn(model_val.A, Wv)
        Wc = opt_rotate(model_cond)
        model_cond.A .= rotate_amn(model_cond.A, Wc)
    end

    if run_maxloc
        @info "Run max localization"
        println()
        model_val.A .= max_localize(model_val)
        model_cond.A .= max_localize(model_cond)
    end

    # Write files
    seedname_val = new_seedname(seedname, outdir_val)
    write_model(seedname_val, model_val)

    seedname_cond = new_seedname(seedname, outdir_cond)
    write_model(seedname_cond, model_cond)

    if rotate_unk
        dir = dirname(seedname)
        outdir_val = dirname(seedname_val)
        outdir_cond = dirname(seedname_cond)
        split_unk(dir, Uv, Uc, outdir_val, outdir_cond)
    end

    return nothing
end

function new_seedname(seedname::String, subdir::String)
    outdir = joinpath(dirname(seedname), subdir)
    !isdir(outdir) && mkdir(outdir)
    return joinpath(outdir, basename(seedname))
end
