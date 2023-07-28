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

- `kpoints`: a length-`n_kpoints` vector of kpoints in fractional coordinates
- `operator_k`: a TB operator on a kgrid, which should be a length-`n_kpoints`
    vector of scalers or matrices or matrices of 3-vector
- `Rspace`: R-space domain for the R-vectors. Since only ``\\mathbf{R}`-vectors
    are used, it can be
    - a [`BareRspace`](@ref)
    - a [`WignerSeitzRspace`](@ref)
    - a [`MDRSRspace`](@ref)

!!! note

    The user should make sure that the `kpoints` are uniformlly distributed,
    since doing a Fourier transform on a non-uniformlly distributed kpoints
    mostly is meaningless in our context, i.e., it will not give us a good
    tight-bindign operator in R-space.
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
    Rspace::AbstractRspace,
)
    nkpts = length(kpoints)
    @assert nkpts > 0 "empty kpoints"
    @assert length(operator_k) == nkpts "operator has wrong n_kpoints"
    nRvecs = n_Rvectors(Rspace)
    @assert length(operator_R) == nRvecs "operator_R has wrong n_Rvectors"

    for (R, Oᴿ) in zip(Rspace.Rvectors, operator_R)
        # clean buffer
        fill!(Oᴿ, zero(eltype(Oᴿ)))
        for (k, Oᵏ) in zip(kpoints, operator_k)
            Oᴿ .+= exp(-im * 2π * dot(k, R)) * Oᵏ
        end
        Oᴿ ./= nkpts
    end
    return nothing
end

function fourier(
    kpoints::AbstractVector, operator_k::AbstractVector, Rspace::AbstractRspace
)
    @assert length(operator_k) > 0 "empty operator_k"
    # force type to complex
    T_op = typeof(complex(first(operator_k[1])))
    size_op = size(operator_k[1])
    operator_R = [zeros(T_op, size_op) for _ in 1:n_Rvectors(Rspace)]
    fourier!(operator_R, kpoints, operator_k, Rspace)
    return operator_R
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

- `Rspace`: the R-space domain, must be a [`BareRspace`](@ref).
- `operator_R`: R-space operator, can be
    - a length-`n_Rvectors` vector of scalers/matrices/matrices of vectors
    - a [`TBOperator`](@ref), then the `Rspace` argument can be omitted
- `kpoints`: kpoints to be interpolated, fractional coordinates, can be nonuniform

!!! note

    The `Rspace` must be a [`BareRspace`](@ref), note also the `Rspace` of
    [`OperatorRspace`](@ref) is also constrained to be [`BareRspace`](@ref).
    This constraint not only simplifies a lot of code, but also makes the
    [`invfourier`](@ref) function much faster.
    One can use [`simplify`](@ref) to convert [`MDRSRspace`](@ref) or
    [`WignerSeitzRspace`](@ref) to [`BareRspace`](@ref) before calling
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
    Rspace::BareRspace,
    operator_R::AbstractVector,
    kpoints::AbstractVector,
)
    nRvecs = n_Rvectors(Rspace)
    @assert length(operator_R) == nRvecs "operator has wrong n_Rvectors"
    nkpts = length(kpoints)
    @assert nkpts > 0 "empty kpoints"
    @assert length(operator_k) == nkpts "operator_k has wrong n_kpoints"

    for (k, Oᵏ) in zip(kpoints, operator_k)
        # clean buffer
        fill!(Oᵏ, zero(eltype(Oᵏ)))
        for (R, Oᴿ) in zip(Rspace.Rvectors, operator_R)
            Oᵏ .+= exp(im * 2π * dot(k, R)) * Oᴿ
        end
    end
    return nothing
end

@inline function invfourier!(
    operator_k::AbstractVector, tb::TBOperator, kpoints::AbstractVector
)
    return invfourier!(operator_k, tb.Rspace, tb.operator, kpoints)
end

function invfourier(Rspace::BareRspace, operator_R::AbstractVector, kpoints::AbstractVector)
    @assert length(operator_R) > 0 "empty operator_R"
    # force type to complex
    T_op = typeof(complex(first(operator_R[1])))
    size_op = size(operator_R[1])
    operator_k = [zeros(T_op, size_op) for _ in 1:length(kpoints)]
    invfourier!(operator_k, Rspace, operator_R, kpoints)
    return operator_k
end

@inline function invfourier(tb::TBOperator, kpoints::AbstractVector)
    return invfourier(tb.Rspace, tb.operator, kpoints)
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
- `kpoints`: kpoints fractional coordinates, should be a uniform grid. Can be
    an iterator with each element a 3-vector, e.g.,
        - a Vector of `Vec3`
- `Rvectors`: fourier frequencies, fractional coordinates, should be iterator
    with each element a 3-vector, e.g.
        - a Vector of `Vec3`
        - a `BareRspace`
        - a `WignerSeitzRspace`
        - a `MDRSRspace`
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
- `Rvectors`: fourier frequencies, fractional coordinates, should be iterator
    with each element a 3-vector, e.g,
    - a `Vector` of `Vec3`
    - a `BareRspace`
- `kpoints`: kpoints to be interpolated, a Vector, each element is a fractional
    coordinates, can be nonuniform
"""
function invfourier(f::Function, Rvectors, kpoints)
    for (ik, k) in enumerate(kpoints)
        for (iR, R) in enumerate(Rvectors)
            phase = exp(im * 2π * (k ⋅ R))
            f(iR, ik, phase)
        end
    end
end
