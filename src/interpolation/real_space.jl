using NearestNeighbors: KDTree, knn

# Internal adapters for the real-space-scheme seam. Both own the quotient-cell
# weights and expose the same three operations to construction.jl:
# `_quotient_vectors`, `_representative_vectors`, and
# `_apply_real_space_selection`. Legacy R-space types are deliberately absent.

struct _WignerSeitzSelection
    quotient_vectors::Vector{Vec3{Int}}
    degeneracies::Vector{Int}

    function _WignerSeitzSelection(
            quotient_vectors::AbstractVector,
            degeneracies::AbstractVector{<:Integer},
        )
        length(quotient_vectors) == length(degeneracies) ||
            throw(DimensionMismatch("R-vector and degeneracy counts differ"))
        all(>(0), degeneracies) ||
            throw(ArgumentError("R-vector degeneracies must be positive"))
        return new(Vector{Vec3{Int}}(quotient_vectors), Vector{Int}(degeneracies))
    end
end

struct _MinimumDistanceContribution
    quotient_index::Int
    row::Int
    column::Int
    denominator::Int
end

struct _MinimumDistanceSelection
    quotient_vectors::Vector{Vec3{Int}}
    representative_vectors::Vector{Vec3{Int}}
    contributions::Vector{Vector{_MinimumDistanceContribution}}
end

_quotient_vectors(selection::_WignerSeitzSelection) = selection.quotient_vectors
_quotient_vectors(selection::_MinimumDistanceSelection) = selection.quotient_vectors
_representative_vectors(selection::_WignerSeitzSelection) = selection.quotient_vectors
_representative_vectors(selection::_MinimumDistanceSelection) =
    selection.representative_vectors

function _validate_real_space_inputs(lattice, grid_size)
    size(lattice) == (3, 3) || throw(ArgumentError("lattice must be a 3 × 3 matrix"))
    length(grid_size) == 3 || throw(ArgumentError("k-point grid size must have length three"))
    all(>(0), grid_size) || throw(ArgumentError("k-point grid sizes must be positive"))
    return nothing
end

function _wigner_seitz_selection(
        lattice::AbstractMatrix,
        grid_size;
        atol::Real,
        max_cell::Integer,
    )
    _validate_real_space_inputs(lattice, grid_size)

    origin = Vec3(0, 0, 0)
    quotient_cell, _ = make_supercell(
        [origin], [0:(size - 1) for size in grid_size]
    )
    candidates, translations = make_supercell(
        quotient_cell, [((-max_cell):max_cell) * size for size in grid_size]
    )
    candidates = sort_points(candidates)
    candidate_cartesian = map(vector -> lattice * vector, candidates)
    translations = unique(translations)
    translation_cartesian = map(vector -> lattice * vector, translations)

    tree = KDTree(translation_cartesian)
    number_neighbors = min(8, length(translation_cartesian))
    nearest_indices, nearest_distances = knn(
        tree, candidate_cartesian, number_neighbors, true
    )
    origin_index = findfirst(==(origin), translations)
    isnothing(origin_index) && error("real-space search does not contain the origin")

    quotient_vectors = Vec3{Int}[]
    degeneracies = Int[]
    for (candidate_index, vector) in enumerate(candidate_cartesian)
        nearest_index = nearest_indices[candidate_index][1]
        nearest_distance = nearest_distances[candidate_index][1]
        if nearest_index != origin_index &&
                abs(nearest_distance - norm(vector)) >= atol
            continue
        end

        push!(quotient_vectors, candidates[candidate_index])
        degeneracy = count(
            distance -> abs(distance - nearest_distance) < atol,
            nearest_distances[candidate_index],
        )
        push!(degeneracies, degeneracy)
    end
    return _WignerSeitzSelection(quotient_vectors, degeneracies)
end

function _minimum_distance_translations(
        lattice::AbstractMatrix,
        grid_size,
        fractional_centers::AbstractVector,
        quotient_vectors::AbstractVector;
        atol::Real,
        max_cell::Integer,
    )
    isempty(fractional_centers) &&
        throw(ArgumentError("minimum-distance selection needs Wannier centers"))
    number_wannier = length(fractional_centers)
    number_quotient_vectors = length(quotient_vectors)

    search_cell = max_cell + 1
    candidates, translations = make_supercell(
        [Vec3(0, 0, 0)],
        [((-search_cell):search_cell) * size for size in grid_size],
    )
    candidate_cartesian = map(vector -> lattice * vector, candidates)
    translations = unique(translations)
    translation_cartesian = map(vector -> lattice * vector, translations)
    tree = KDTree(translation_cartesian)
    number_neighbors = min(8, length(translation_cartesian))
    origin_index = findfirst(==(Vec3(0, 0, 0)), translations)
    isnothing(origin_index) && error("minimum-distance search does not contain the origin")

    selected_translations = Array{Vector{Vec3{Int}}, 3}(
        undef, number_wannier, number_wannier, number_quotient_vectors
    )
    for (quotient_index, quotient_vector) in enumerate(quotient_vectors)
        for column in 1:number_wannier, row in 1:number_wannier
            displacement =
                fractional_centers[column] + quotient_vector - fractional_centers[row]
            displaced_candidates = map(
                vector -> vector + lattice * displacement, candidate_cartesian
            )
            nearest_indices, nearest_distances = knn(
                tree, displaced_candidates, number_neighbors, true
            )

            selected_indices = Int[]
            for candidate_index in eachindex(displaced_candidates)
                indices = nearest_indices[candidate_index]
                distances = nearest_distances[candidate_index]
                nearest_index = indices[1]
                if nearest_index != origin_index
                    origin_neighbor = findfirst(==(origin_index), indices)
                    isnothing(origin_neighbor) && continue
                    abs(distances[1] - distances[origin_neighbor]) < atol || continue
                end
                push!(selected_indices, candidate_index)
            end
            isempty(selected_indices) && error(
                "minimum-distance search found no representative for " *
                    "Wannier pair ($row, $column) and quotient vector $quotient_vector",
            )
            selected_translations[row, column, quotient_index] =
                candidates[selected_indices]
        end
    end
    return selected_translations
end

function _minimum_distance_selection(
        wigner_seitz::_WignerSeitzSelection,
        selected_translations::AbstractArray{<:AbstractVector, 3},
        translation_degeneracies::AbstractArray{<:Integer, 3},
    )
    quotient_vectors = wigner_seitz.quotient_vectors
    number_quotient_vectors = length(quotient_vectors)
    size(selected_translations, 3) == number_quotient_vectors ||
        throw(DimensionMismatch("translation and quotient-vector counts differ"))
    number_wannier = size(selected_translations, 1)
    size(selected_translations, 2) == number_wannier ||
        throw(DimensionMismatch("minimum-distance translation matrices must be square"))
    size(translation_degeneracies) == size(selected_translations) ||
        throw(DimensionMismatch("translation and degeneracy shapes differ"))

    representative_vectors = Vec3{Int}[]
    representative_indices = Dict{Vec3{Int}, Int}()
    contributions = Vector{Vector{_MinimumDistanceContribution}}()
    for (quotient_index, quotient_vector) in enumerate(quotient_vectors)
        quotient_degeneracy = wigner_seitz.degeneracies[quotient_index]
        for column in 1:number_wannier, row in 1:number_wannier
            translations = selected_translations[row, column, quotient_index]
            translation_degeneracy =
                translation_degeneracies[row, column, quotient_index]
            translation_degeneracy > 0 ||
                throw(ArgumentError("minimum-distance translation sets cannot be empty"))
            length(translations) == translation_degeneracy || throw(
                ArgumentError(
                    "the number of minimum-distance translations must match its degeneracy",
                ),
            )
            denominator = quotient_degeneracy * translation_degeneracy
            for translation in translations
                vector = Vec3{Int}(quotient_vector + translation)
                representative_index = get!(representative_indices, vector) do
                    push!(representative_vectors, vector)
                    push!(contributions, _MinimumDistanceContribution[])
                    return length(representative_vectors)
                end
                push!(
                    contributions[representative_index],
                    _MinimumDistanceContribution(
                        quotient_index, row, column, denominator
                    ),
                )
            end
        end
    end
    return _MinimumDistanceSelection(
        copy(quotient_vectors), representative_vectors, contributions
    )
end

function _real_space_selection(
        scheme::WignerSeitz,
        lattice,
        grid_size,
        fractional_centers;
        atol = scheme.atol,
        max_cell = scheme.max_cell,
    )
    return _wigner_seitz_selection(lattice, grid_size; atol, max_cell)
end

function _real_space_selection(
        scheme::MinimumDistance,
        lattice,
        grid_size,
        fractional_centers;
        atol = scheme.atol,
        max_cell = scheme.max_cell,
    )
    wigner_seitz = _wigner_seitz_selection(lattice, grid_size; atol, max_cell)
    translations = _minimum_distance_translations(
        lattice,
        grid_size,
        fractional_centers,
        wigner_seitz.quotient_vectors;
        atol,
        max_cell,
    )
    translation_degeneracies = map(length, translations)
    return _minimum_distance_selection(
        wigner_seitz, translations, translation_degeneracies
    )
end

function _apply_real_space_selection(
        selection::_WignerSeitzSelection,
        quotient_coefficients::AbstractArray,
    )
    number_dimensions = ndims(quotient_coefficients)
    size(quotient_coefficients, number_dimensions) == length(selection.quotient_vectors) ||
        throw(DimensionMismatch("coefficient and quotient-vector counts differ"))
    selected = similar(quotient_coefficients)
    for quotient_index in eachindex(selection.quotient_vectors)
        source = selectdim(quotient_coefficients, number_dimensions, quotient_index)
        destination = selectdim(selected, number_dimensions, quotient_index)
        destination .= source ./ selection.degeneracies[quotient_index]
    end
    return selected
end

function _apply_real_space_selection(
        selection::_MinimumDistanceSelection,
        quotient_coefficients::AbstractArray,
    )
    number_dimensions = ndims(quotient_coefficients)
    number_wannier = size(quotient_coefficients, 1)
    size(quotient_coefficients, 2) == number_wannier ||
        throw(DimensionMismatch("the two Wannier matrix dimensions must be equal"))
    size(quotient_coefficients, number_dimensions) == length(selection.quotient_vectors) ||
        throw(DimensionMismatch("coefficient and quotient-vector counts differ"))

    components = Tuple(size(quotient_coefficients)[3:(end - 1)])
    number_components = prod(components; init = 1)
    packed_source = reshape(
        quotient_coefficients,
        number_wannier,
        number_wannier,
        number_components,
        length(selection.quotient_vectors),
    )
    packed_selected = zeros(
        eltype(quotient_coefficients),
        number_wannier,
        number_wannier,
        number_components,
        length(selection.representative_vectors),
    )
    for (representative_index, entries) in enumerate(selection.contributions)
        for entry in entries
            source = view(
                packed_source,
                entry.row,
                entry.column,
                :,
                entry.quotient_index,
            )
            destination = view(
                packed_selected,
                entry.row,
                entry.column,
                :,
                representative_index,
            )
            destination .+= source ./ entry.denominator
        end
    end
    return reshape(
        packed_selected,
        number_wannier,
        number_wannier,
        components...,
        length(selection.representative_vectors),
    )
end
