"""
A submodule to allow easy loading of artifacts as `Model`s.
"""
module Datasets

# need this for lazy artifacts
using LazyArtifacts

using Artifacts

# for pretty printing directory tree
using FileTrees

using ..Wannier: read_w90

export list_datasets, load_dataset, show_dataset, @dataset_str

function list_datasets()
    artifacts_toml = Artifacts.find_artifacts_toml(@__DIR__)
    artifact_dict = Artifacts.load_artifacts_toml(artifacts_toml)
    return collect(keys(artifact_dict))
end

function load_dataset(name)
    # the directory of the artifact
    dir = Artifacts.@artifact_str(name)

    # I assume the prefix (`prefix.mmn`) is the same as the folder name
    prefix = name

    return read_w90(joinpath(dir, prefix))
end

function show_dataset(name)
    return FileTree(Artifacts.@artifact_str(name))
end

"""
Works like `Artifacts.@artifact_str`, but the `Artifacts.toml` points to the
`Wannier.jl/Artifacts.toml`.

Defined as a new macro to avoid overwriting the original `@artifact_str`.
Cannot directly `export @artifact_str` because it can not locate `Artifacts.toml`.

# Example

```julia
using Wannier.Datasets
dataset"Si2"
```
"""
macro dataset_str(s)
    @eval Artifacts.@artifact_str $s
end

end
