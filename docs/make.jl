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
