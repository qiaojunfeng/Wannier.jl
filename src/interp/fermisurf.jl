using ProgressBars: ProgressBar

@doc raw"""
    fermi_surface(Rvectors::RV, H::AbstractArray{Complex{T},3}; n_k::KT) where {T<:Real, RV<:Union{RVectors{T},RVectorsMDRS{T}},KT<:Union{AbstractVector{Int},Integer}}

Interpolate Fermi surface.

# Arguments
- `Rvectors`: `RVectors` or `RVectorsMDRS`
- `H`: `n_wann * n_wann * n_rvecs`, Hamiltonian in R space
- `n_k`: integer or 3-vector, number of interpolated kpoints along three directions

# Return
- `kpoints`: `3 * n_kx * n_ky * n_kz`, interpolated kpoints in fractional coordinates
- `E`: `n_wann * n_kx * n_ky * n_kz`, interpolated eigenvalues

!!! note

    The output `n_kx = n_k + 1`, since the bxsf format requires general grid, i.e.,
    the last kpoint is the periodic image of the first one.
    This also restores the behavior of `Wannier90`.

!!! note

    For MDRS interpolation, the `H` should be defined on the ``\bm{R}`` vectors
    instead of the MDRSv2 ``\tilde{\bm{R}}`` vectors; the expansion of ``H(\bm{R})``
    to ``H(\tilde{\bm{R}})`` is done internally. This means that if you read the
    `seedname_tb.dat` file, then you can directly pass the `H` to this function.
"""
function fermi_surface(
    Rvectors::RV, H::AbstractArray{Complex{T},3}; n_k::KT
) where {
    T<:Real,RV<:Union{RVectors{T},RVectorsMDRS{T}},KT<:Union{AbstractVector{Int},Integer}
}
    n_wann, _, n_rvecs = size(H)
    n_rvecs == Rvectors.n_rvecs || error("n_rvecs of H != Rvectors.n_rvecs")

    n_kx, n_ky, n_kz = _expand_nk(n_k)
    # kpoints are in fractional coordinates
    kpoints = get_kpoints([n_kx, n_ky, n_kz]; endpoint=true)
    n_kpts = n_kx * n_ky * n_kz

    E = zeros(T, n_wann, n_kpts)

    # By default I use the MDRSv2 for MDRS interpolation, so I expand it here
    if Rvectors isa RVectorsMDRS
        H = mdrs_v1tov2(Rvectors, H)
    end

    println("n_threads: ", Threads.nthreads())

    # Actually doing this once is faster, but there's no multi-threading
    Hwork = zeros(Complex{T}, n_wann, n_wann, n_kpts)
    invfourier!(Hwork, Rvectors, H, kpoints)
    # There are lots of FFTs, so I add a progress bar, although it slows down a bit
    for ik in ProgressBar(1:n_kpts)
        Hᵏ = @view Hwork[:, :, ik]
        # check Hermiticity
        @assert norm(Hᵏ - Hᵏ') < 1e-10
        # diagonalize
        ϵ = eigen(Hᵏ).values
        E[:, ik] = real.(ϵ)
    end

    # # preallocate buffers
    # Hwork = [zeros(Complex{T}, n_wann, n_wann, 1) for _ in 1:Threads.nthreads()]
    # # There are lots of FFTs, so I add a progress bar, although it slows down a bit
    # Threads.@threads for ik in ProgressBar(1:n_kpts)
    #     k = kpoints[:, ik:ik]  # use range in last dimension to keep it 2D
    #     # Hwork[ik] .= invfourier(Rvectors, H, k)
    #     invfourier!(Hwork[Threads.threadid()], Rvectors, H, k)
    #     Hᵏ = @view Hwork[Threads.threadid()][:, :, 1]
    #     # check Hermiticity
    #     @assert norm(Hᵏ - Hᵏ') < 1e-10
    #     # diagonalize
    #     ϵ = eigen(Hᵏ).values
    #     E[:, ik] = real.(ϵ)
    # end

    # The kz increase the fastest in kpoints, reshape them to (n_kx, n_ky, n_kz)
    kpoints = reshape(kpoints, 3, n_kz, n_ky, n_kx)
    E = reshape(E, n_wann, n_kz, n_ky, n_kx)
    # and permutedims
    kpoints = permutedims(kpoints, (1, 4, 3, 2))
    E = permutedims(E, (1, 4, 3, 2))

    return kpoints, E
end

"""
    fermi_surface(model; n_k)

Interpolate Fermi surface.

# Arguments
- `model`: [`InterpModel`](@ref)

See also [`fermi_surface`](@ref fermi_surface(Rvectors::RV, H::AbstractArray{Complex{T},3}; n_k::KT)).
"""
function fermi_surface(model::InterpModel; n_k)
    return fermi_surface(model.kRvectors.Rvectors, model.H; n_k=n_k)
end

function _expand_nk(n_k::T) where {T<:Union{AbstractVector{Int},Integer}}
    length(n_k) in [1, 3] || error("n_k must be an integer or a vector of length 3")

    # bxsf need general grid, i.e., the last kpoint is the periodic image of the first one
    # so I increase the n_k by 1
    if length(n_k) == 3
        n_kx, n_ky, n_kz = [n + 1 for n in n_k]
    else
        n_kx = n_ky = n_kz = n_k + 1
    end

    return n_kx, n_ky, n_kz
end
