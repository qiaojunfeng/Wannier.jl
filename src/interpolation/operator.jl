export simplify

"""
An asbtract type representing a tight-binding operator.

Since julia does not permit inheritance for fields, we need to ask the developers
to implement the following fields when introducing a new type of tight-binding
operator, see [`HamiltonianRspace`](@ref) for an example.

# Fields
- `domain`: the domain on which the operator is defined, can be
    - a [`BareRspaceDomain`](@ref) for R-space operators
    - a [`AbstractKpointContainer`](@ref) for k-space operators
- `operator`: the actual operator defined on the `domain`.
    - for R-space, it is a length-`n_Rvectors` vector,
    - for k-space, it is a length-`n_kpoints` vector.
    Each element type `T` can be
    - a scaler
    - a matrix of scaler, e.g., `Matrix{ComplexF64}` for Hamiltonian
    - a matrix of vector, e.g., `Matrix{Vec3{ComplexF64}}` for position
        operator, spin operator, etc.
"""
abstract type AbstractTBOperator end

function n_wannier(tb_operator::AbstractTBOperator)
    return isempty(tb_operator.operator) ? 0 : size(tb_operator.operator[1], 1)
end

# the domain is so simple so we expose it to the user
function Base.propertynames(tb::AbstractTBOperator)
    return Tuple([collect(fieldnames(typeof(tb.domain))); collect(fieldnames(typeof(tb)))])
end

function Base.getproperty(tb::AbstractTBOperator, sym::Symbol)
    type_domain = typeof(getfield(tb, :domain))
    if sym ∈ fieldnames(type_domain)
        return getfield(getfield(tb, :domain), sym)
    else
        # fallback
        return getfield(tb, sym)
    end
end

function Base.show(io::IO, ::MIME"text/plain", tb::AbstractTBOperator)
    @printf(io, "TB operator type    :  %s\n", nameof(typeof(tb)))
    show(io, MIME"text/plain"(), tb.domain)
    println(io)
    @printf(io, "n_wannier   =  %d", n_wannier(tb))
end

"""Index using `i` of domain"""
function Base.getindex(tb::AbstractTBOperator, i::Integer)
    return tb.operator[i]
end
function Base.getindex(tb::AbstractTBOperator, r::UnitRange{<:Integer})
    return tb.operator[r]
end

Base.lastindex(tb::AbstractTBOperator) = lastindex(tb.operator)
Base.length(tb::AbstractTBOperator) = length(tb.operator)

function Base.zeros(tb::AbstractTBOperator)
    return typeof(tb)(tb.domain, [zeros(tb.operator[1]) for _ in 1:n_Rvectors(tb)])
end
function Base.fill!(tb::AbstractTBOperator, x)
    for O in tb.operator
        if eltype(O) <: AbstractArray
            for Oi in O
                fill!(Oi, x)
            end
        else
            O .= x
        end
    end
end

"""
An asbtract type representing a tight-binding operator in R-space.

See [`AbstractTBOperator`](@ref) for the fields that must be implemented.

In this case, the `domain` field refers to the R-space domain (or called
R-vectors), and its type must be [`BareRspaceDomain`](@ref) since
[`invfourier`](@ref) only accepts that type.
"""
abstract type AbstractOperatorRspace <: AbstractTBOperator end

n_Rvectors(tb::AbstractOperatorRspace) = n_Rvectors(tb.domain)
real_lattice(tb::AbstractOperatorRspace) = real_lattice(tb.domain)

"""Index using 3-vector of R-vectors coordinates"""
function Base.getindex(tb::AbstractOperatorRspace, R::Vec3{Int})
    return getindex(tb, R[1], R[2], R[3])
end

"""Index using 3-vector of R-vectors coordinates"""
function Base.getindex(tb::AbstractOperatorRspace, x::Integer, y::Integer, z::Integer)
    iR = tb.domain.xyz_iR[x, y, z]
    if iszero(iR)
        throw(BoundsError(tb, [x, y, z]))
    end
    return tb.operator[iR]
end

"""
    $(SIGNATURES)

Simplify the R-vectors of a tight-binding operator, such that the inverse
Fourier transform is a simple sum (and is faster).

- If the R-space domain is a `MDRSRspaceDomain`, expand the R-vectors to a new
    set of R-vectors `R̃ = R + T`, where `T` is the translation vectors;
    also divide the operator by the degeneracies of `R` and `T` vectors.
- If the R-space domain is a `WSRspaceDomain`, divide the operator by the
    R-vector degeneracy.
- If the R-space domain is a `BareRspaceDomain`, do nothing.

# Return
- a `BareRspaceDomain` with simplified Rvectors
- a vector for the simplified `operator` defined on the `BareRspaceDomain`
"""
function simplify(Rdomain::AbstractRspaceDomain, operator::AbstractVector) end

function simplify(Rdomain::MDRSRspaceDomain, operator::AbstractVector)
    # expanded R-vectors
    bare_Rvectors = Vector{Vec3{Int}}()
    # simplified operator by absorbing R and T degeneracies
    bare_operator = similar(operator, 0)
    op_type = eltype(operator[1])
    op_size = size(operator[1])
    # generate expanded R vectors, which contains all the R+T
    for iR in 1:n_Rvectors(Rdomain)
        R = Rdomain.Rvectors[iR]
        Tvecs = Rdomain.Tvectors[iR]
        Nᴿ = Rdomain.n_Rdegens[iR]
        for n in axes(Tvecs, 2)
            for m in axes(Tvecs, 1)
                Nᵀ = Rdomain.n_Tdegens[iR][m, n]
                for iT in 1:Nᵀ
                    RT = R .+ Tvecs[m, n][iT]
                    i = findfirst(x -> x == RT, bare_Rvectors)
                    if isnothing(i)
                        push!(bare_Rvectors, RT)
                        Oᴿᵀ = zeros(op_type, op_size)
                        Oᴿᵀ[m, n] = operator[iR][m, n] / (Nᴿ * Nᵀ)
                        push!(bare_operator, Oᴿᵀ)
                    else
                        bare_operator[i][m, n] += operator[iR][m, n] / (Nᴿ * Nᵀ)
                    end
                end
            end
        end
    end
    bare_Rdomain = BareRspaceDomain(Rdomain.lattice, bare_Rvectors)
    return bare_Rdomain, bare_operator
end

function simplify(Rdomain::WSRspaceDomain, operator::AbstractVector)
    # absorb R degeneracies into operator
    bare_operator = map(zip(operator, Rdomain.n_Rdegens)) do (O, Nᴿ)
        O ./ Nᴿ
    end
    bare_Rdomain = BareRspaceDomain(Rdomain.lattice, Rdomain.Rvectors)
    return bare_Rdomain, bare_operator
end

function simplify(Rdomain::BareRspaceDomain, operator::AbstractVector)
    return Rdomain, operator
end

"""
An asbtract type representing a tight-binding operator in k-space.

# Fields
Must implement the following fields, see [`HamiltonianKspace`](@ref) for an example:
- `domain`: the k-space domain on which the operator is defined, can be either one of
    - a [`KpointList`](@ref) for non-uniformly distributed kpoints
    - a [`KpointGrid`](@ref) for uniformlly distributed kpoints
- `operator`: the tight-binding operator defined on the `klist`

!!! note

    We differentiate betwenn `KpointList` and `KpointGrid` because we limit
    [`fourier`](@ref) only on `KpointGrid` type, but [`invfourier`] can
    accept both types. This is because doing a Fourier transform on a
    non-uniformlly distributed kpoints mostly is meaningless in our context,
    i.e., it will not give us a good TB operator in R-space.
"""
abstract type AbstractOperatorKspace{K<:AbstractKpointContainer} <: AbstractTBOperator end

n_kpoints(tb::AbstractOperatorKspace) = n_kpoints(tb.domain)
reciprocal_lattice(tb::AbstractOperatorKspace) = reciprocal_lattice(tb.domain)

"""
    $(SIGNATURES)

Sort a Rspace operator (e.g. Hamiltonian ``H(\\mathbf{R})``) by the norm of
``\\mathbf{R}``-vectors.

Since the operator decays ~expoentially with ``\\mathbf{R}``, the resulting
operator should be in descending order.

# Return
- `idx`: the index to sort the operator, i.e. `Oᴿ[idx]` is the sorted operator
- `normR`: the norm of R vectors, `normR[idx]` is the sorted `R` vectors in ascending order
- `normO`: the norm of the operator, `normO[idx]` is the sorted `Oᴿ` in descending order

# Example
```julia
hamiltonian, position = read_w90_tb("mos2")
idx, normR, normH = Wannier.sort_by_Rnorm(hamiltonian)
```
"""
function sort_by_Rnorm(operator_R::AbstractOperatorRspace)
    normO = map(norm, operator_R.operator)
    normR = map(operator_R.Rvectors) do R
        norm(operator_R.lattice * R)
    end
    idx = sortperm(normR)

    @printf("# No.    ||R|| (Å)    ||H(R)||\n")
    for i in 1:20
        @printf("%2d       %.5f       %.5f\n", i - 1, normR[idx[i]], normO[idx[i]])
    end
    return idx, normR, normO
end

"""
    $(SIGNATURES)

Cut R-space operator (e.g. Hamiltonian) by the given cutoff radius.

# Arguments
- `Rcut`: cutoff radius in angstrom, the operator `O(R)` with R-vector norm
    larger than `Rcut` will removed.

# Example
```julia
hamiltonian, position = read_w90_tb("mos2")
H_cut = Wannier.cut(hamiltonian)
```
"""
function cut(operator_R::AbstractOperatorRspace, Rcut::Real)
    println("Sorting by R norm...")
    _, normR, _ = sort_by_Rnorm(operator_R)
    println("")
    idx_keep = normR .<= Rcut
    println(typeof(operator_R.domain))
    Rdomain = typeof(operator_R.domain)(
        operator_R.domain.lattice, deepcopy(operator_R.Rvectors[idx_keep])
    )
    return typeof(operator_R)(Rdomain, deepcopy(operator_R.operator[idx_keep]))
end
