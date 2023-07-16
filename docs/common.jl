# Copied from
# https://github.com/thchr/Brillouin.jl/blob/fad88c5b6965fe4bd59e725ea60655348d36ce0f/docs/make.jl#L4
# ---------------------------------------------------------------------------------------- #
# make PlotlyJS plots showable in ```@example ``` blocks, following the approach suggested
# in https://github.com/fredrikekre/Literate.jl/issues/126
using PlotlyJS

struct HTMLPlot
    p
    h::Int # desired display height in pixels
end

HTMLPlot(p) = HTMLPlot(p, 400)

# the Documenter.jl build folder for storing the final HTML pages
const ROOT_DIR = joinpath(@__DIR__, "build")
# the folder for saving HTML plots, needs to be inside the Documenter.jl build folder
const PLOT_DIR = joinpath(ROOT_DIR, "plots")

function Base.show(io::IO, ::MIME"text/html", p::HTMLPlot)
    mkpath(PLOT_DIR)
    path = joinpath(PLOT_DIR, string(hash(p) % UInt32, ".html"))
    PlotlyJS.savefig(p.p, path; format="html")
    return print(
        io,
        "<object type=\"text/html\" data=\"/$(relpath(path, ROOT_DIR))\" style=\"width:100%;height:$(p.h)px;\"></object>",
    )
end
# ---------------------------------------------------------------------------------------- #

# add nbviewer and mybinder links
# copied from DFTK
# https://github.com/JuliaMolSim/DFTK.jl/blob/master/docs/make.jl
function add_badges(str)
    badges = [
        "[![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/examples/@__NAME__.ipynb)",
        "[![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/examples/@__NAME__.ipynb)",
    ]

    # Find the Header and insert the badges right below
    splitted = split(str, "\n")
    idx = findfirst(startswith.(splitted, "# # "))
    idx === nothing && error("Example files must start with # #")
    insert!(splitted, idx + 1, "#md # " * badges[1])
    insert!(splitted, idx + 2, "#md # " * badges[2])
    insert!(splitted, idx + 3, "#md #")
    return join(splitted, "\n")
end
