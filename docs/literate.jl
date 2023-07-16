using Literate

include("common.jl")

println("\n**** Process examples with Literate.jl ****")

# dir of the original literate jl script
LITERATE_SRCDIR = joinpath(@__DIR__, "literate")
# dir of literate generated markdown/ipynb files
LITERATE_OUTDIR = joinpath(@__DIR__, "src")

for (root, _, files) in walkdir(LITERATE_SRCDIR), file in files
    splitext(file)[2] == ".jl" || continue

    println("\n==== Processing $file ====")

    srcpath = joinpath(root, file)
    outdir = splitdir(replace(srcpath, LITERATE_SRCDIR => LITERATE_OUTDIR))[1]

    # generate markdown which will be executed by Documenter.jl
    Literate.markdown(srcpath, outdir)

    # I skip the execution of the notebook, because
    # 1. it increases the build time
    # 2. somehow ipynb does not show the plots correctly, e.g. bands, WFs, etc.
    # 3. random numbers during execution might cause the notebook output to be different
    # 4. I will let the user download an empty notebook, so that at least they will run once :-)
    Literate.notebook(srcpath, outdir; execute=false)

    # This generates a cleansed version w/o comments
    # Literate.script(srcpath, outdir)
end

println("\n**** Literate.jl finished ****\n")
