using Documenter
using Literate
using Wannier

TUTORIALS_SRCDIR = joinpath(@__DIR__, "src/tutorials_src")
TUTORIALS_OUTDIR = joinpath(@__DIR__, "src/tutorials")
for md in readdir(TUTORIALS_SRCDIR)
    file = joinpath(TUTORIALS_SRCDIR, md)
    Literate.markdown(file, TUTORIALS_OUTDIR)
    # Literate.notebook(file, TUTORIALS_OUTDIR)
    Literate.script(file, TUTORIALS_OUTDIR)
end

makedocs(;
    sitename="Wannier.jl",
    authors="Junfeng Qiao and contributors.",
    modules=[Wannier],
    # the `example` blocks in the tutorials need correct path to read files
    workdir=joinpath(@__DIR__, "../tutorials/tutorials/"),
    pages=[
        "Home" => "index.md",
        "Quick start" => "quickstart.md",
        # the tutorials will be processed by Literate
        "Tutorial" => ["Maximal localization" => "tutorials/1-maxloc.md"],
        "Theory" => ["Normalization" => "theory/normalization.md"],
        "API" => [
            "Utilities" => "api/util.md",
            "Input/Output" => "api/io.md",
            "B vectors" => "api/bvector.md",
            "Model" => "api/model.md",
            "Wannierize" => "api/wannierize.md",
            "Interpolation" => "api/interpolation.md",
            "Real space" => "api/realspace.md",
            "Command line" => "api/cli.md",
        ],
    ],
)
