using Literate

# Process tutorials by Literate.jl

for dir in readdir(TUTORIALS_SRCDIR)
    # skip folders starting with underscore
    # Some times I put work-in-progress folders there but do not want them to be
    # processed by Literate.jl and Documenter.jl.
    isdir(joinpath(TUTORIALS_SRCDIR, dir)) && !startswith(dir, "_") || continue

    println("* Processing tutorial: $dir *")

    # quick hack to skip generated folders during development
    # startswith(dir, [Char('0' + i) for i in 1:9]) && continue

    # jl = filter(x -> endswith(x, ".jl"), readdir(joinpath(TUTORIALS_SRCDIR, dir)))
    # I assume the tutorial file is named as `tutorial.jl`
    jl = "tutorial.jl"

    file = joinpath(TUTORIALS_SRCDIR, dir, jl)
    isfile(file) || error("tutorial file not found: $file")

    outdir = joinpath(TUTORIALS_OUTDIR, dir)

    # generate markdown which will be executed by Documenter.jl
    Literate.markdown(file, outdir)

    # the Documenter.jl execution needs input files, so I copy them to
    # the build folder where the markdown will be executed
    for f in readdir(joinpath(TUTORIALS_SRCDIR, dir))
        # skip the tutorial file
        f == jl && continue

        src = joinpath(TUTORIALS_SRCDIR, dir, f)
        dstdir = joinpath(TUTORIALS_BUILDDIR, dir)
        dst = joinpath(dstdir, f)

        mkpath(dstdir)
        cp(src, dst; force=true, follow_symlinks=true)
    end

    # I skip the execution of the notebook, because
    # 1. it increases the build time
    # 2. somehow ipynb does not show the plots correctly, e.g. bands, WFs, etc.
    # 3. random numbers during execution might cause the notebook output to be different
    # 4. I will let the user download an empty notebook, so that at least they will run once :-)
    Literate.notebook(file, outdir; execute=false)

    # generate a cleansed version w/o comments, so that user can run in CLI
    Literate.script(file, outdir)
end

# process the WannierPlots examples

for jl in readdir(PLOTS_EXAMPLES_SRCDIR)
    # only process jl scripts
    isfile(joinpath(PLOTS_EXAMPLES_SRCDIR, jl)) && endswith(jl, ".jl") || continue

    file = joinpath(PLOTS_EXAMPLES_SRCDIR, jl)
    println("* Processing WannierPlots example: $jl *")

    # quick hack to skip generated folders during development
    # startswith(jl, [Char('0' + i) for i in 1:3]) && continue

    # save the markdown just inside the srcdir
    Literate.markdown(file, PLOTS_EXAMPLES_SRCDIR)

    # I skip the execution of the notebook
    Literate.notebook(file, PLOTS_EXAMPLES_SRCDIR; execute=false)
end

# Copy the input files to the build folder, so that Documenter.jl can correctly
# execute the markdown files.

# since the whole build dir is deleted, I need to create it again
mkpath(PLOTS_EXAMPLES_BUILDDIR)

for dir in ["3-band", "4-realspace", "8-fermi_surface"]
    src = joinpath(TUTORIALS_SRCDIR, dir)
    dst = joinpath(PLOTS_EXAMPLES_BUILDDIR, dir)
    cp(src, dst; force=true)
end

println("* Literate.jl finished *\n")
