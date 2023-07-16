using Literate

include("common.jl")

println("\n**** Process examples with Literate.jl ****")

# dir of the original literate jl script
LITERATE_SRCDIR = joinpath(@__DIR__, "literate")
# dir of literate generated markdown/ipynb files
LITERATE_OUTDIR = joinpath(@__DIR__, "src")

PROCESS_ALL_SCRIPTS = false
if PROCESS_ALL_SCRIPTS
    # this process all the jl scripts
    src_files = []
    for (root, _, files) in walkdir(LITERATE_SRCDIR), file in files
        append!(src_files, joinpath(root, file))
    end
else
    # while this only process the files in `EXAMPLES`, so that I can comment out
    # some examples in `make.jl` to speed up the build
    src_files = map(EXAMPLES) do (_, file)
        file = joinpath(LITERATE_SRCDIR, replace(file, ".md" => ".jl"))
    end
end

for file in src_files
    splitext(file)[2] == ".jl" || continue

    println("\n==== Processing $file ====")

    outdir = splitdir(replace(file, LITERATE_SRCDIR => LITERATE_OUTDIR))[1]

    # generate markdown which will be executed by Documenter.jl
    Literate.markdown(file, outdir; preprocess=add_badges)

    # I skip the execution of the notebook, because
    # 1. it increases the build time
    # 2. somehow ipynb does not show the plots correctly, e.g. bands, WFs, etc.
    # 3. random numbers during execution might cause the notebook output to be different
    # 4. I will let the user download an empty notebook, so that at least they will run once :-)
    Literate.notebook(file, outdir; execute=false)

    # This generates a cleansed version w/o comments
    # Literate.script(file, outdir)
end

println("\n**** Literate.jl finished ****\n")
