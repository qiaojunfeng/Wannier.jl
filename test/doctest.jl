using Documenter
using Wannier

DocMeta.setdocmeta!(Wannier, :DocTestSetup, :(using Wannier); recursive = true)
doctest(
    Wannier
    # fix=true,  # update all the output in `jldoctest`
)
