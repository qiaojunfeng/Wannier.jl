using Documenter
using Literate
using Wannier

TUTORIALS_SRCDIR = joinpath(@__DIR__, "src/tutorials_src")
TUTORIALS_OUTDIR = joinpath(@__DIR__, "src/tutorials")

# Copied from https://github.com/thchr/Brillouin.jl/blob/fad88c5b6965fe4bd59e725ea60655348d36ce0f/docs/make.jl#L4
# ---------------------------------------------------------------------------------------- #
# make PlotlyJS plots showable in ```@example ``` blocks, following the approach suggested
# in https://github.com/fredrikekre/Literate.jl/issues/126
using PlotlyJS
struct HTMLPlot
    p
    h::Int # desired display height in pixels
end
HTMLPlot(p) = HTMLPlot(p, 400)
const ROOT_DIR = joinpath(@__DIR__, "build/tutorials")
const PLOT_DIR = joinpath(ROOT_DIR, "plots")
function Base.show(io::IO, ::MIME"text/html", p::HTMLPlot)
    mkpath(PLOT_DIR)
    path = joinpath(PLOT_DIR, string(hash(p) % UInt32, ".html"))
    PlotlyJS.savefig(p.p, path; format="html")
    return print(
        io,
        "<object type=\"text/html\" data=\"../$(relpath(path, ROOT_DIR))\" style=\"width:100%;height:$(p.h)px;\"></object>",
    )
end
# ---------------------------------------------------------------------------------------- #

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
        "Tutorial" => [
            "Maximal localization" => "tutorials/1-maxloc.md",
            "Disentanglement" => "tutorials/2-disentangle.md",
            "Band structure" => "tutorials/3-band.md",
        ],
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
