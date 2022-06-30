#!/usr/bin/env julia
using Wannier


function main()
    # seedname = "silicon"
    seedname = "/home/jqiao/git/Wannier.jl/test/fixtures/silicon"

    # Input AMN is Silicon s,p projection
    model = read_seedname(seedname)

    read_chk = true
    if read_chk
        # Get max localized gauge from chk file
        chk = read_chk("$seedname.chk.fmt")
        # We replace the initial projection by the "good" max loc gauge
        model.A = get_amn(chk)
    else
        # You can also use disentangle to get a good gauge from initial projection
        model.A = disentangle(model)
    end

    n_val = model.n_wann ÷ 2

    # UNK files for plotting WFs
    rotate_unk = false
    splitted = split_wannierize(model, n_val, rotate_unk)
    if rotate_unk
        model_val, model_cond, Uv, Uc = splitted
    else
        model_val, model_cond = splitted
    end

    # Write files
    outdir_val = "val"
    seedname_val = new_seedname(seedname, outdir_val)
    write_model(seedname_val, model_val)

    outdir_cond = "cond"
    seedname_cond = new_seedname(seedname, outdir_cond)
    write_model(seedname_cond, model_cond)

    if rotate_unk
        dir = dirname(seedname)
        outdir_val = dirname(seedname_val)
        outdir_cond = dirname(seedname_cond)
        split_unk(dir, Uv, Uc, outdir_val, outdir_cond)
    end

    nothing
end


function new_seedname(seedname::String, subdir::String)
    outdir = joinpath(dirname(seedname), subdir)
    !isdir(outdir) && mkdir(outdir)
    joinpath(outdir, basename(seedname))
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
