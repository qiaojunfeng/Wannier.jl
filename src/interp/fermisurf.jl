using ProgressBars: ProgressBar

@doc raw"""
    fermi_surface(Rvectors, H; n_k)

Interpolate Fermi surface.

# Arguments
- `Rvectors`: `RVectors` or `RVectorsMDRS`
- `H`: `n_wann * n_wann * n_r̃vecs`, Hamiltonian in R space
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
function fermi_surface(hami::TBHamiltonian; n_k::KT
) where {
    KT<:Union{AbstractVector{Int},Integer}
}
    nwann = n_wann(hami)
    n_rvecs = length(hami)
    
    n_kx, n_ky, n_kz = _expand_nk(n_k)
    @printf("Interpolation grid: %d %d %d\n", n_kx, n_ky, n_kz)

    # kpoints are in fractional coordinates
    kpoints = get_kpoints([n_kx, n_ky, n_kz]; endpoint=true)
    n_kpts = n_kx * n_ky * n_kz

    println("n_threads: ", Threads.nthreads())

    #     E[:, ik] = real.(ϵ)
    # end
    hk = HamiltonianKGrid(hami, kpoints)

    # TODO Don't forget to redo the correct order in bsxf
    # The kz increase the fastest in kpoints, reshape them to (n_kx, n_ky, n_kz)
    # kpoints = reshape(kpoints, 3, n_kz, n_ky, n_kx)
    # E = reshape(E, nwann, n_kz, n_ky, n_kx)
    # # and permutedims
    # kpoints = permutedims(kpoints, (1, 4, 3, 2))
    # E = permutedims(E, (1, 4, 3, 2))

    return kpoints, hk.eigvals
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
