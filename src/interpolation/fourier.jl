export fourier!, invfourier!, fourier, invfourier

"""
    $(SIGNATURES)

Fourier transform an operator from k-space to R-space.

```math
O_{mn}(\\mathbf{R}) = \\frac{1}{N_{\\mathbf{k}}}
\\sum_{\\mathbf{k}} \\exp(-i {\\mathbf{k}} \\mathbf{R}) O_{mn}(\\mathbf{k}),
```
where
- ``N_{\\mathbf{k}}``: the total number of kpoints

# Arguments
There are multiple methods for this function, the following are the explanations
of all possible arguments, but not all of them are needed for each method:

- `Rdomain`: R-space domain for the R-vectors. Since only ``\\mathbf{R}`-vectors
    are used, it can be
    - a [`BareRspaceDomain`](@ref)
    - a [`MDRSRspaceDomain`](@ref)
    - a [`WSRspaceDomain`](@ref)
- `kpoints`: a length-`n_kpoints` vector of kpoints in fractional coordinates
- `operator_k`: a TB operator on a kgrid, which can be
    - a length-`n_kpoints` vector of scalers/matrices/matrices of vectors
    - a [`AbstractOperatorKspace`](@ref), but must be defined on a
        [`KpointGrid`](@ref) so that the Fourier transform is on a
        uniformlly-distributed kpoint grid.

!!! note

    The user should make sure that the `kpoints` are uniformlly distributed,
    otherwise the forward Fourier transform does not make sense.
"""
function fourier end

"""
In-place version of [`fourier`](@ref).

# Arguments
- `operator_R`: the output R-space operator, will be overwritten

See [`fourier`](@ref) for other arguments.
"""
function fourier! end

function fourier!(
    operator_R::AbstractVector,
    kpoints::AbstractVector,
    operator_k::AbstractVector,
    Rdomain::AbstractRspaceDomain,
)
    nkpts = length(kpoints)
    @assert nkpts > 0 "empty kpoints"
    @assert length(operator_k) == nkpts "operator has wrong n_kpoints"
    nRvecs = n_Rvectors(Rdomain)
    @assert length(operator_R) == nRvecs "operator_R has wrong n_Rvectors"

    for (R, Oᴿ) in zip(Rdomain.Rvectors, operator_R)
        Oᴿ .= 0  # clean buffer
        for (k, Oᵏ) in zip(kpoints, operator_k)
            Oᴿ .+= exp(-im * 2π * dot(k, R)) * Oᵏ
        end
        Oᴿ ./= nkpts
    end
    return nothing
end

@inline function fourier!(
    operator_R::AbstractVector,
    operator_k::AbstractOperatorKspace{K},
    Rdomain::AbstractRspaceDomain,
) where {K<:KpointGrid}
    return fourier!(operator_R, operator_k.domain.kpoints, operator_k.operator, Rdomain)
end

function fourier(
    kpoints::AbstractVector, operator_k::AbstractVector, Rdomain::AbstractRspaceDomain
)
    @assert length(operator_k) > 0 "empty operator_k"
    T_op = complex(eltype(operator_k[1]))
    size_op = size(operator_k[1])
    operator_R = [zeros(T_op, size_op) for _ in 1:n_Rvectors(Rdomain)]
    fourier!(operator_R, kpoints, operator_k, Rdomain)
    return operator_R
end

@inline function fourier(
    operator_k::AbstractOperatorKspace{K}, Rdomain::AbstractRspaceDomain
) where {K<:KpointGrid}
    return fourier(operator_k.domain.kpoints, operator_k.operator, Rdomain)
end

"""
    $(SIGNATURES)

Inverse Fourier transform operator from R-space to k-space.

```math
O_{mn}(\\mathbf{k}) = \\sum_{\\mathbf{R}}
\\exp(i \\mathbf{k} \\mathbf{R}) O_{mn}(\\mathbf{R}).
```

# Arguments
There are multiple methods for this function, the following are the explanations
of all possible arguments, but not all of them are needed for each method:

- `kpoints`: kpoints to be interpolated, fractional coordinates, can be nonuniform
- `Rdomain`: the R-space domain, must be a [`BareRspaceDomain`](@ref).
- `operator_R`: R-space operator, can be
    - a length-`n_Rvectors` vector of scalers/matrices/matrices of vectors
    - a [`AbstractOperatorRspace`](@ref)

!!! note

    The `Rdomain` must be a [`BareRspaceDomain`](@ref), note also the
    `Rdomain` of [`AbstractOperatorRspace`](@ref) is also constrained to be
    [`BareRspaceDomain`](@ref).
    This constraint not only simplifies a lot of code, but also makes the
    [`invfourier`](@ref) function much faster.
    One can use [`simplify`](@ref) to convert [`MDRSRspaceDomain`](@ref) or
    [`WSRspaceDomain`](@ref) to [`BareRspaceDomain`](@ref) before calling
    this function.
"""
function invfourier end

"""
In-place version of [`invfourier`](@ref).

# Arguments
- `operator_k`: the output k-space operator, will be overwritten

See [`invfourier`](@ref) for other arguments.
"""
function invfourier! end

function invfourier!(
    operator_k::AbstractVector,
    Rdomain::BareRspaceDomain,
    operator_R::AbstractVector,
    kpoints::AbstractVector,
)
    nRvecs = n_Rvectors(Rdomain)
    @assert length(operator_R) == nRvecs "operator has wrong n_Rvectors"
    nkpts = length(kpoints)
    @assert nkpts > 0 "empty kpoints"
    @assert length(operator_k) == nkpts "operator_k has wrong n_kpoints"

    for (k, Oᵏ) in zip(kpoints, operator_k)
        Oᵏ .= 0  # clean buffer
        for (R, Oᴿ) in zip(Rdomain.Rvectors, operator_R)
            Oᵏ .+= exp(im * 2π * dot(k, R)) * Oᴿ
        end
    end
    return nothing
end

@inline function invfourier!(
    operator_k::AbstractVector, operator_R::AbstractOperatorRspace, kpoints::AbstractVector
)
    return invfourier!(operator_k, operator_R.domain, operator_R.operator, kpoints)
end

@inline function invfourier!(
    operator_k::AbstractVector,
    operator_R::AbstractOperatorRspace,
    kdomain::AbstractKpointContainer,
)
    return invfourier!(operator_k, operator_R.domain, operator_R.operator, kdomain.kpoints)
end

@inline function invfourier!(
    operator_k::AbstractOperatorKspace, operator_R::AbstractOperatorRspace
)
    return invfourier!(operator_k.operator, operator_R, operator_k.domain)
end

function invfourier(
    Rdomain::BareRspaceDomain, operator_R::AbstractVector, kpoints::AbstractVector
)
    @assert length(operator_R) > 0 "empty operator_R"
    T_op = complex(eltype(operator_R[1]))
    size_op = size(operator_R[1])
    operator_k = [zeros(T_op, size_op) for _ in 1:length(kpoints)]
    invfourier!(operator_k, Rdomain, operator_R, kpoints)
    return operator_k
end

@inline function invfourier(operator_R::AbstractOperatorRspace, kpoints::AbstractVector)
    return invfourier(operator_R.domain, operator_R.operator, kpoints)
end

@inline function invfourier(
    operator_R::AbstractOperatorRspace, kdomain::AbstractKpointContainer
)
    return invfourier(operator_R, kdomain.kpoints)
end

@inline function invfourier(
    operator_R::AbstractOperatorRspace, operator_k::AbstractOperatorKspace
)
    return invfourier(operator_R, operator_k.domain)
end

"""
    $(SIGNATURES)

Fourier transform from the k-space to the R-space.

# Arguments
- `f`: function to be called, should have signature `f(ik, iR, phase)`:
    - `ik`: index of kpoint
    - `iR`: index of Rvector
    - `phase`: the phase factor, which is
        ``\\exp(-2\\pi i \\mathbf{k} \\mathbf{R}) / N_{\\mathbf{k}}``
        where ``N_{\\mathbf{k}}`` is the number of kpoints
- `kpoints`: kpoints to be interpolated, fractional coordinates, should be a
    uniform grid. Should be an iterator with each element a 3-vector, e.g.,
        - a Vector of `Vec3`
        - [`AbstractKpointContainer`](@ref)
- `Rvectors`: fourier frequencies, fractional coordinates, can be nonuniform.
    Should be iterator with each element a 3-vector.
"""
function fourier(f::Function, kpoints, Rvectors)
    nkpts = length(kpoints)
    for (iR, R) in enumerate(Rvectors)
        for (ik, k) in enumerate(kpoints)
            phase = exp(-im * 2π * (k ⋅ R)) / nkpts
            f(ik, iR, phase)
        end
    end
end

"""
    $(SIGNATURES)

Inverse Fourier transform from R-space to k-space.

# Arguments
- `f`: function to be called, should have signature `f(iR, ik, phase)`:
    - `ik`: index of kpoint
    - `iR`: index of Rvector
    - `phase`: the phase factor, which is
        ``\\exp(2\\pi i \\mathbf{k} \\mathbf{R})``
- `Rvectors`: fourier frequencies, fractional coordinates, can be nonuniform.
    Should be iterator with each element a 3-vector.
- `kpoints`: kpoints to be interpolated, fractional coordinates, should be a
    uniform grid. Should be an iterator with each element a 3-vector, e.g.,
        - a Vector of `Vec3`
        - [`AbstractKpointContainer`](@ref)
"""
function invfourier(f::Function, Rvectors, kpoints)
    for (ik, k) in enumerate(kpoints)
        for (iR, R) in enumerate(Rvectors)
            phase = exp(im * 2π * (k ⋅ R))
            f(iR, ik, phase)
        end
    end
end
