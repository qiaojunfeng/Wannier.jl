export TBOperator, simplify

"""
A struct representing a tight-binding operator in R-space.

# Fields
$(FIELDS)

!!! note

    We limit the type of `Rspace` to only [`BareRspace`](@ref), because
    this simplifies a lot of codes and compute faster.
    To covert a [`WignerSeitzRspace`](@ref) or [`MDRSRspace`](@ref)
    to [`BareRspace`](@ref), see [`simplify`](@ref).
    Note also that [`invfourier`](@ref) only accepts [`BareRspace`](@ref).
"""
struct TBOperator{M<:AbstractMatrix}
    """a concise name for the operator"""
    name::String

    """the R-space domain (or called R-vectors) on which the operator is defined"""
    Rspace::BareRspace

    """the actual operator defined on the `Rspace`, should be a
    length-`n_Rvectors` vector, where each element type `T` can be
    - a scaler
    - a matrix of scaler, e.g., `Matrix{ComplexF64}` for Hamiltonian
    - a matrix of vector, e.g., `Matrix{MVec3{ComplexF64}}` for position
        operator, spin operator, etc."""
    operator::Vector{M}
end

function n_wannier(tb::TBOperator)
    return isempty(tb.operator) ? 0 : size(tb.operator[1], 1)
end
n_Rvectors(tb::TBOperator) = n_Rvectors(tb.Rspace)
real_lattice(tb::TBOperator) = real_lattice(tb.Rspace)

# the Rspace is so simple so we expose it to the user
function Base.propertynames(tb::TBOperator)
    return Tuple([
        collect(fieldnames(typeof(tb.Rspace)))
        collect(fieldnames(typeof(tb)))
    ])
end

function Base.getproperty(tb::TBOperator, sym::Symbol)
    type_R = typeof(getfield(tb, :Rspace))
    if sym ∈ fieldnames(type_R)
        return getfield(getfield(tb, :Rspace), sym)
    else
        # fallback
        return getfield(tb, sym)
    end
end

function Base.show(io::IO, ::MIME"text/plain", tb::TBOperator)
    @printf(io, "Tight-binding operator name  :  %s\n", tb.name)
    show(io, MIME"text/plain"(), tb.Rspace)
    println(io)
    @printf(io, "n_wannier   =  %d", n_wannier(tb))
end

"""Index using `i` of Rspace"""
function Base.getindex(tb::TBOperator, i::Integer)
    return tb.operator[i]
end
function Base.getindex(tb::TBOperator, r::UnitRange{<:Integer})
    return tb.operator[r]
end

Base.lastindex(tb::TBOperator) = lastindex(tb.operator)
Base.length(tb::TBOperator) = length(tb.operator)

function Base.iterate(tb::TBOperator, state=1)
    if state > length(tb.operator)
        return nothing
    else
        return (tb.operator[state], state + 1)
    end
end

"""Index using 3-vector of R-vectors coordinates"""
function Base.getindex(tb::TBOperator, R::Vec3{Int})
    return getindex(tb, R[1], R[2], R[3])
end

"""Index using 3-vector of R-vectors coordinates"""
function Base.getindex(tb::TBOperator, x::Integer, y::Integer, z::Integer)
    iR = tb.Rspace.xyz_iR[x, y, z]
    if iszero(iR)
        throw(BoundsError(tb, [x, y, z]))
    end
    return tb.operator[iR]
end

function Base.zeros(tb::TBOperator)
    @assert length(tb) > 0 "empty operator"
    T = eltype(tb.operator[1])
    s = size(tb.operator[1])
    return TBOperator("zeros", tb.Rspace, [zeros(T, s) for _ in 1:length(tb)])
end

function Base.fill!(tb::TBOperator, x)
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

function Base.isapprox(a::TBOperator, b::TBOperator; kwargs...)
    return isapprox_struct(a, b; kwargs...)
end

"""
    $(SIGNATURES)

Simplify the R-vectors of a tight-binding operator, such that the inverse
Fourier transform is a simple sum (and is faster).

- If the R-space domain is a `MDRSRspace`, expand the R-vectors to a new
    set of R-vectors `R̃ = R + T`, where `T` is the translation vectors;
    also divide the operator by the degeneracies of `R` and `T` vectors.
- If the R-space domain is a `WignerSeitzRspace`, divide the operator by the
    R-vector degeneracy.
- If the R-space domain is a `BareRspace`, do nothing.

# Return
- a `BareRspace` with simplified Rvectors
- a vector for the simplified `operator` defined on the `BareRspace`
"""
function simplify end

function simplify(Rspace::MDRSRspace, operator::AbstractVector)
    # expanded R-vectors
    bare_Rvectors = Vector{Vec3{Int}}()
    # simplified operator by absorbing R and T degeneracies
    bare_operator = similar(operator, 0)
    op_type = eltype(operator[1])
    op_size = size(operator[1])
    # generate expanded R vectors, which contains all the R+T
    for iR in 1:n_Rvectors(Rspace)
        R = Rspace.Rvectors[iR]
        Tvecs = Rspace.Tvectors[iR]
        Nᴿ = Rspace.n_Rdegens[iR]
        for n in axes(Tvecs, 2)
            for m in axes(Tvecs, 1)
                Nᵀ = Rspace.n_Tdegens[iR][m, n]
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
    bare_Rspace = BareRspace(Rspace.lattice, bare_Rvectors)
    return bare_Rspace, bare_operator
end

function simplify(Rspace::WignerSeitzRspace, operator::AbstractVector)
    # absorb R degeneracies into operator
    bare_operator = map(zip(operator, Rspace.n_Rdegens)) do (O, Nᴿ)
        O ./ Nᴿ
    end
    bare_Rspace = BareRspace(Rspace.lattice, Rspace.Rvectors)
    return bare_Rspace, bare_operator
end

function simplify(Rspace::BareRspace, operator::AbstractVector)
    return Rspace, operator
end

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
function sort_by_Rnorm(tb::TBOperator)
    normO = map(norm, tb.operator)
    normR = map(tb.Rvectors) do R
        norm(tb.lattice * R)
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
function cut(tb::TBOperator, Rcut::Real)
    println("Sorting by R norm...")
    _, normR, _ = sort_by_Rnorm(tb)
    println("")
    idx_keep = normR .<= Rcut
    Rspace = BareRspace(tb.Rspace.lattice, deepcopy(tb.Rvectors[idx_keep]))
    return TBOperator(tb.name, Rspace, deepcopy(tb.operator[idx_keep]))
end

"""An abstract type that interpolate physical quantities from some
[`TBOperator`](@ref)s

Each subtype should define a function
`function (interp::name_of_Interpolator)(kpoints::AbstractVector{<:AbstractVector})`
that returns the interpolated physical quantities for the given list of `kpoints`.

Since it interpolates back to Bloch gauge, almost always the 1st field is
`hamiltonian::TBHamiltonian`.
"""
abstract type AbstractTBInterpolator <: Function end

#=
I cannot define a function like this, because this will cause method ambiguity:
    - the KPathInterpolant is also a AbstractVector{<:AbstractVector}
    - each concrete interpolator type defines a function for AbstractVector{<:AbstractVector}

@inline function (interp::AbstractTBInterpolator)(kpi::KPathInterpolant; kwargs...)
    kpoints = get_kpoints(kpi)
    return interp(kpoints; kwargs...)
end
=#

@inline function (interp::AbstractTBInterpolator)(
    kpoint::AbstractVector{<:Real}, args...; kwargs...
)
    # unwrap results
    result = interp([kpoint], args...; kwargs...)
    if length(result) == 1
        return result[1]
    else
        return Tuple(map(x -> x[1], result))
    end
end

n_wannier(interp::AbstractTBInterpolator) = n_wannier(interp.hamiltonian)
n_Rvectors(interp::AbstractTBInterpolator) = n_Rvectors(interp.hamiltonian)
real_lattice(interp::AbstractTBInterpolator) = real_lattice(interp.hamiltonian)

function Base.show(io::IO, ::MIME"text/plain", interp::AbstractTBInterpolator)
    @printf(io, "Tight-binding interpolator name  :  %s\n", nameof(interp))
    ops = []
    itps = []
    for p in propertynames(interp)
        f = getproperty(interp, p)
        if f isa TBOperator
            push!(ops, f.name)
        elseif f isa AbstractTBInterpolator
            push!(itps, nameof(typeof(f)))
        else
            continue
        end
    end
    if length(itps) > 0
        println(io, "\nList of contained interpolators:")
        for it in itps
            println(io, "  ", it)
        end
    end
    if length(ops) > 0
        println(io, "\nList of contained operators:")
        for op in ops
            println(io, "  ", op)
        end
    end
    println(io, "\nSummary:")
    length(itps) > 0 && @printf(io, "  n_interpolators =  %d\n", length(itps))
    length(ops) > 0 && @printf(io, "  n_operators     =  %d\n", length(ops))
    @printf(io, "  n_Rvectors      =  %d\n", n_Rvectors(interp))
    @printf(io, "  n_wannier       =  %d", n_wannier(interp))
end
