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

const ROOT_DIR = joinpath(@__DIR__, "build")
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
