using Literate

include("common.jl")

println("\n**** Process examples with Literate.jl ****")

# Requires variable `EXAMPLES` in caller script
for (_, md) in EXAMPLES
    println("\n==== Processing $md ====")

    endswith(md, ".md") || continue
    # replace
    jl = replace(md, r"\.md$" => ".jl")

    file = joinpath(@__DIR__, "src", jl)
    isfile(file) || error("tutorial file not found: $file")

    outdir = dirname(file)

    # generate markdown which will be executed by Documenter.jl
    Literate.markdown(file, outdir)

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
