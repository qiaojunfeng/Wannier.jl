"""
A submodule to allow easy loading of artifacts as `Model`s.
"""
module Datasets

# need this for lazy artifacts
using LazyArtifacts

# do not bring exported stuff into scope since I need to extend the @artifact_str macro
using Artifacts: Artifacts

# for pretty printing directory tree
using FileTrees

using ..Wannier: read_w90

export list_datasets, load_dataset, show_dataset, @artifact_str

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
Load artifact from `Wannier.jl/Artifacts.toml`.

Cannot directly `export @artifact_str` because it can not locate `Artifacts.toml`.
"""
macro artifact_str(s)
    @eval Artifacts.@artifact_str $s
end

end
