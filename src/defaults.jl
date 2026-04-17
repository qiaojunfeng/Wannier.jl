using Dates: now
using Printf: @sprintf

default_header() = @sprintf "Created by Wannier.jl %s" string(now())

"""Default is to use all the bands."""
default_exclude_bands() = Vector{Int}()

"""
Wannier90's default value for its input parameter `kmesh_tol`.
To reproduce the same results when generating ``\\mathbf{b}``-vectors.
"""
default_w90_kmesh_tol() = 1.0e-6

"""
Wannier90's default value for its input parameter `search_shells`.
To reproduce the same results when generating ``\\mathbf{b}``-vectors.
"""
default_w90_bvectors_search_shells() = 36

"""
Wannier90's internal constant for checking parallel bvectors.
See [`are_parallel`](@ref).
"""
default_w90_bvectors_check_parallel_atol() = 1.0e-6

"""
Wannier90's internal constant for checking the sigular value when computing
bvector bweights. See [`compute_bweights`](@ref).
"""
default_w90_bvectors_singular_value_atol() = 1.0e-5

"""
Wannier90's internal constant for sorting bvector order.
See [`sort_bvectors`](@ref).
"""
default_w90_bvectors_sort_supercell_atol() = 1.0e-8

"""
    default_entangled_bands(nkpoints, nbands)

Default entangled bands, i.e., all participating in disentanglement.
"""
function default_entangled_bands end

function default_entangled_bands(nkpoints::Integer, nbands::Integer)
    return trues(nbands, nkpoints)
end

function default_entangled_bands(gauges::AbstractArray{<:Any, 3})
    return trues(size(gauges, 1), size(gauges, 3))
end

"""Wannier90's default number of kpoints along the 1st segment of kpath."""
default_w90_kpath_num_points() = 100

"""Wannier90's default tolerance for searching R-vectors."""
default_w90_ws_distance_tol() = 1.0e-5

"""Wannier90's default search range for R-vectors."""
default_w90_ws_search_size() = 3

"""Wannier90's default threshold for degenerate eigenvalues in berry-related calculations."""
default_w90_berry_degen_tol() = 1.0e-4

"""Wannier90's default switch for using degenerate perturbation theory."""
default_w90_berry_use_degen_pert() = false

"""Wannier90's default switch for enforcing Hermiticity of position operator."""
default_w90_berry_position_force_hermiticity() = true

"""Wannier90's default switch for enforcing Hermiticity of
position*hamiltonian*position operator, for orbital magnetization."""
default_w90_berry_duHdu_force_hermiticity() = true
