export BandStructure, band_structure

"""Band energies and plotting coordinates along a labeled reciprocal-space path."""
struct BandStructure{T <: Real, K <: AbstractVector, E <: AbstractMatrix}
    kpoints::K
    path_coordinate::Vector{T}
    symmetry_point_indices::Vector{Int}
    symmetry_point_labels::Vector{String}
    band_energy::E
end

Base.length(bands::BandStructure) = length(bands.kpoints)

"""
    band_structure(model, kpath)

Interpolate the band energy along a `KPath`, retaining its linear path
coordinate and high-symmetry labels.
"""
function band_structure(model::InterpolationModel, kpath::KPath)
    path_coordinate, symmetry_point_indices, symmetry_point_labels = linear_path(kpath)
    result = interpolate(model, kpath.points, BandEnergy())
    return BandStructure(
        copy(kpath.points),
        path_coordinate,
        symmetry_point_indices,
        symmetry_point_labels,
        result.band_energy,
    )
end
