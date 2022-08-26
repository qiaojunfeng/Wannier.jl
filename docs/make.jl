using Documenter
using Wannier

makedocs(;
    sitename="Wannier.jl",
    authors="Junfeng Qiao and contributors.",
    modules=[Wannier],
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Theory" => ["Normalization" => "theory/normalization.md"],
        "API" => "api.md",
    ],
)
