using LinearAlgebra

@doc raw"""
    fourier(kRvectors::KRVectors{T,RV}, Oᵏ::Array{Complex{T},3})
    where {T<:Real,RV<:RVectors{T}}

Fourier transform operator from k space to R space for Wiger-Seitz interpolation.

```math
X_{mn}(\bm{R}) = \frac{1}{N_{\bm{k}}}
\sum_{\bm{k}} \exp(-i {\bm{k}} \bm{R}) X_{mn}(\bm{k}),
```
where `N_{\bm{k}}` is the total number of kpoints.

# Arguments
- `kRvectors`: for WS interpolation
- `Oᵏ`: `n_wann * n_wann * n_kpts`, operator matrix in k space
"""
function fourier(
    kRvectors::KRVectors{T,RV}, Oᵏ::Array{Complex{T},3}
) where {T<:Real,RV<:RVectors{T}}
    n_wann = size(Oᵏ, 1)
    n_kpts = kRvectors.n_kpts
    n_rvecs = kRvectors.n_rvecs
    @assert size(Oᵏ, 1) == size(Oᵏ, 2)
    @assert size(Oᵏ, 3) == n_kpts

    Oᴿ = zeros(Complex{T}, n_wann, n_wann, n_rvecs)

    for ir in 1:n_rvecs
        R = kRvectors.Rvectors.R[:, ir]
        for ik in 1:n_kpts
            k = kRvectors.kpoints[:, ik]
            fac = exp(-im * 2π * dot(k, R))
            Oᴿ[:, :, ir] += fac * Oᵏ[:, :, ik]
        end
    end

    Oᴿ ./= n_kpts
    return Oᴿ
end

@doc raw"""
    invfourier(kRvectors::KRVectors{T,RV}, Oᴿ::Array{Complex{T},3}, kpoints::AbstractMatrix{T})
    where {T<:Real,RV<:RVectors{T}}

Inverse Fourier transform operator from R space to k space for WS interpolation.

```math
X_{mn}(\bm{k}) = \sum_{\bm{R}} \frac{1}{N_{\bm{R}}}
\exp(i \bm{k} \bm{R}) X_{mn}(\bm{R}),
```
where ``N_{\bm{R}}`` is the degeneracy of R vectors (not the total number of R vectors).

# Arguments
- `Oᴿ`: `n_wann * n_wann * n_rvecs`, real space operator matrix
- `kpoints`: `3 * n_kpts`, kpoints to be interpolated, fractional coordinates, can be nonuniform
"""
function invfourier(
    kRvectors::KRVectors{T,RV}, Oᴿ::Array{Complex{T},3}, kpoints::AbstractMatrix{T}
) where {T<:Real,RV<:RVectors{T}}
    n_wann = size(Oᴿ, 1)
    n_kpts = size(kpoints, 2)
    n_rvecs = kRvectors.n_rvecs
    @assert size(Oᴿ, 1) == size(Oᴿ, 2)
    @assert size(Oᴿ, 3) == n_rvecs

    Oᵏ = zeros(Complex{T}, n_wann, n_wann, n_kpts)

    for ik in 1:n_kpts
        k = kpoints[:, ik]
        for ir in 1:n_rvecs
            R = kRvectors.Rvectors.R[:, ir]
            fac = exp(im * 2π * dot(k, R))
            fac /= kRvectors.Rvectors.N[ir]
            Oᵏ[:, :, ik] += fac * Oᴿ[:, :, ir]
        end
    end
    return Oᵏ
end

@doc raw"""
    _fourier_mdrs_v1(kRvectors::KRVectors{T,RV}, Oᵏ::Array{Complex{T},3})
    where {T<:Real,RV<:RVectorsMDRS{T}}

Fourier transform operator from k space to R space for MDRSv1.

The same as the forward Fourier transform of WS interpolation,
```math
X_{mn}(\bm{R}) = \frac{1}{N_{\bm{k}}}
\sum_{\bm{k}} \exp(-i \bm{k} \bm{R}) X_{mn}(\bm{k}),
```
where `N_{\bm{k}}` is the total number of kpoints.

See also [`_invfourier_mdrs_v1`](@ref _invfourier_mdrs_v1).

# Arguments
- `kRvectors`: for MDRS interpolation
- `Oᵏ`: `n_wann * n_wann * n_kpts`, operator matrix in k space
"""
function _fourier_mdrs_v1(
    kRvectors::KRVectors{T,RV}, Oᵏ::Array{Complex{T},3}
) where {T<:Real,RV<:RVectorsMDRS{T}}
    # this is the same as Wigner-Seitz fourier
    kR = KRVectors{T,RVectors{T}}(
        kRvectors.lattice,
        kRvectors.kgrid,
        kRvectors.kpoints,
        kRvectors.k_xyz,
        kRvectors.xyz_k,
        kRvectors.Rvectors.Rvectors,
        kRvectors.recip_lattice,
        kRvectors.n_kpts,
    )
    Oᴿ = fourier(kR, Oᵏ)
    return Oᴿ
end

@doc raw"""
    _fourier_mdrs_v2(kRvectors::KRVectors{T,RV}, Oᵏ::Array{Complex{T},3})
    where {T<:Real,RV<:RVectorsMDRS{T}}

Fourier transform operator from k space to R space for MDRSv2.

To differentiate from the ``{X}_{mn}(\bm{R})`` of MDRSv1, I use
``\widetilde{X}_{mn}(\widetilde{ \bm{R} })`` for MDRSv2,
```math
\widetilde{X}_{mn}(\widetilde{ \bm{R} }) =
\sum_{ \bm{R} } \frac{1}{ N_{ \bm{R} } \mathcal{N}_{ mn \bm{R} } } X_{mn}(\bm{R})
\sum_{j=1}^{\mathcal{N}_{mn \bm{R} }}
\delta_{ \widetilde{ \bm{R} }, \bm{R} + \bm{T}_{ mn \bm{R} }^{(j)} },
```
where ``N_{\bm{R}}`` is the degeneracy of R vectors (not the total number of R vectors),
and ``\mathcal{N}_{ mn \bm{R} }`` is the degeneracy of ``\bm{T}_{ mn \bm{R} }`` vectors.

See also [`_invfourier_mdrs_v2`](@ref _invfourier_mdrs_v2).

# Arguments
- `kRvectors`: for MDRS interpolation
- `Oᵏ`: `n_wann * n_wann * n_kpts`, operator matrix in k space

!!! note

    This is a faster version of MDRS, by removing some for loops.
    But need to expand R vectors to `R + T`, so different from the R vectors
    in `Wannier90` output `seeedname_tb.dat`.
"""
function _fourier_mdrs_v2(
    kRvectors::KRVectors{T,RV}, Oᵏ::Array{Complex{T},3}
) where {T<:Real,RV<:RVectorsMDRS{T}}
    # first generate O(R) where R is just the WS interpolation R vectors
    Oᴿ = _fourier_mdrs_v1(kRvectors, Oᵏ)
    # now expand O(R) to O(R̃) where R̃ is the MDRS expanded R vectors
    n_r̃vecs = kRvectors.n_r̃vecs
    n_wann = size(Oᵏ, 1)
    Oᴿ̃ = zeros(Complex{T}, n_wann, n_wann, n_r̃vecs)

    for ir̃ in 1:n_r̃vecs
        for (ir, m, n, it) in kRvectors.Rvectors.R̃_RT[ir̃]
            fac = kRvectors.Rvectors.Nᵀ[m, n, ir]
            # I divide here the degeneracy or R vector,
            # so no need to divide again in invfourier
            fac *= kRvectors.Rvectors.N[ir]
            Oᴿ̃[m, n, ir̃] += Oᴿ[m, n, ir] / fac
        end
    end
    return Oᴿ̃
end

"""
    fourier(kRvectors::KRVectors{T,RV}, Oᵏ::Array{Complex{T},3}; version::Symbol=:v2)
    where {T<:Real,RV<:RVectorsMDRS{T}}

Wrapper function for Fourier transform for both MDRS v1 and v2.

# Arguments
- `kRvectors`: for WS interpolation
- `Oᵏ`: `n_wann * n_wann * n_kpts`, operator matrix in k space

# Keyword Arguments
- `version`: `:v1` or `:v2`, default is `:v2`

See also [`_fourier_mdrs_v1`](@ref _fourier_mdrs_v1) and [`_fourier_mdrs_v2`](@ref _fourier_mdrs_v2).
"""
function fourier(
    kRvectors::KRVectors{T,RV}, Oᵏ::Array{Complex{T},3}; version::Symbol=:v2
) where {T<:Real,RV<:RVectorsMDRS{T}}
    version ∈ [:v1, :v2] || error("version must be v1 or v2")
    if version == :v1
        return _fourier_mdrs_v1(kRvectors, Oᵏ)
    else
        return _fourier_mdrs_v2(kRvectors, Oᵏ)
    end
end

@doc raw"""
    _invfourier_mdrs_v1(
        kRvectors::KRVectors{FT,RV}, Oᴿ::Array{Complex{FT},3}, kpoints::AbstractMatrix{FT}
    ) where {FT<:Real,RV<:RVectorsMDRS{FT}}

Inverse Fourier transform operator from R space to k space for MDRSv1 interpolation.

```math
X_{mn}(\bm{k}) = \sum_{\bm{R}}
\frac{1}{ \mathcal{N}_{mn \bm{R}} } X_{mn}(\bm{R})
\sum_{j=1}^{\mathcal{N}_{ mn \bm{R} }}
\exp\left( i \bm{k} \cdot \left( \bm{R} + \bm{T}_{ mn \bm{R} }^{(j)} \right) \right)
```
where ``N_{\bm{R}}`` is the degeneracy of R vectors (not the total number of R vectors),
and ``\mathcal{N}_{ mn \bm{R} }`` is the degeneracy of ``\bm{T}_{ mn \bm{R} }`` vectors.

See also [`_fourier_mdrs_v1`](@ref _fourier_mdrs_v1).

# Arguments
- `Oᴿ`: `n_wann * n_wann * n_rvecs`, real space operator matrix
- `kpoints`: `3 * n_kpts`, kpoints to be interpolated, fractional coordinates, can be nonuniform

!!! note

    This is the MDRSv1 using plain loops.
    It is slower, but reproduce the `Wannier90` behavior,
    and generate the same `seedname_tb.dat` file as `Wannier90`.
    The MDRSv2 has expanded set of R vectors, thus cannot be compared with `tb.dat` file.
"""
function _invfourier_mdrs_v1(
    kRvectors::KRVectors{FT,RV}, Oᴿ::Array{Complex{FT},3}, kpoints::AbstractMatrix{FT}
) where {FT<:Real,RV<:RVectorsMDRS{FT}}
    n_wann = size(Oᴿ, 1)
    n_kpts = size(kpoints, 2)
    n_rvecs = kRvectors.n_rvecs
    @assert size(Oᴿ, 1) == size(Oᴿ, 2)
    @assert size(Oᴿ, 3) == n_rvecs

    Oᵏ = zeros(Complex{FT}, n_wann, n_wann, n_kpts)

    for ik in 1:n_kpts
        k = kRvectors.kpoints[:, ik]
        for ir in 1:n_rvecs
            R = kRvectors.Rvectors.R[:, ir]
            N = kRvectors.Rvectors.N[ir]
            for n in 1:n_wann
                for m in 1:n_wann
                    Nᵀ = kRvectors.Rvectors.Nᵀ[m, n, ir]
                    for T in eachcol(kRvectors.Rvectors.T[m, n, ir])
                        fac = exp(im * 2π * dot(k, (R + T)))
                        fac /= N * Nᵀ
                        Oᵏ[m, n, ik] += fac * Oᴿ[m, n, ir]
                    end
                end
            end
        end
    end
    return Oᵏ
end

@doc raw"""
    _invfourier_mdrs_v2(
        kRvectors::KRVectors{T,RV}, Oᴿ̃::Array{Complex{T},3}, kpoints::AbstractMatrix{T}
    ) where {T<:Real,RV<:RVectorsMDRS{T}}

Inverse Fourier transform operator from R space to k space for MDRSv2 interpolation.

```math
X_{mn}(\bm{k}) = \sum_{ \widetilde{\bm{R}} }
\exp\left( i \bm{k} \widetilde{ \bm{R} } \right) \widetilde{X}_{mn}(\widetilde{\bm{R}})
```

See also [`_fourier_mdrs_v2`](@ref _fourier_mdrs_v2).

# Arguments
- `Oᴿ̃`: `n_wann * n_wann * n_rvecs`, real space operator matrix
- `kpoints`: `3 * n_kpts`, kpoints to be interpolated, fractional coordinates, can be nonuniform
"""
function _invfourier_mdrs_v2(
    kRvectors::KRVectors{T,RV}, Oᴿ̃::Array{Complex{T},3}, kpoints::AbstractMatrix{T}
) where {T<:Real,RV<:RVectorsMDRS{T}}
    n_wann = size(Oᴿ̃, 1)
    n_kpts = size(kpoints, 2)
    n_r̃vecs = kRvectors.Rvectors.n_r̃vecs
    @assert size(Oᴿ̃, 1) == size(Oᴿ̃, 2)
    @assert size(Oᴿ̃, 3) == n_r̃vecs

    Oᵏ = zeros(Complex{T}, n_wann, n_wann, n_kpts)

    # almost the same as WS invfourier, except that
    # the degeneracies are handled during fourier transform
    for ik in 1:n_kpts
        k = kpoints[:, ik]
        for ir̃ in 1:n_r̃vecs
            R̃ = kRvectors.Rvectors.R̃vectors.R[:, ir̃]
            fac = exp(im * 2π * dot(k, R̃))
            Oᵏ[:, :, ik] += fac * Oᴿ̃[:, :, ir̃]
        end
    end
    return Oᵏ
end

"""
    invfourier(
        kRvectors::KRVectors{T,RV}, Oᴿ::Array{Complex{T},3}, kpoints::AbstractMatrix{T};
        version::Symbol=:v2,
    ) where {T<:Real,RV<:RVectorsMDRS{T}}

Wrapper function for inverse Fourier transform for both MDRS v1 and v2.

# Arguments
- `kRvectors`: for MDRS interpolation
- `Oᴿ`: `n_wann * n_wann * n_rvecs`, operator matrix in R space
- `kpoints`: `3 * n_kpts`, kpoints to be interpolated, fractional coordinates, can be nonuniform

# Keyword Arguments
- `version`: `:v1` or `:v2`, default is `:v2`

See also [`_invfourier_mdrs_v1`](@ref _invfourier_mdrs_v1) and
[`_invfourier_mdrs_v2`](@ref _invfourier_mdrs_v2).
"""
function invfourier(
    kRvectors::KRVectors{T,RV},
    Oᴿ::Array{Complex{T},3},
    kpoints::AbstractMatrix{T};
    version::Symbol=:v2,
) where {T<:Real,RV<:RVectorsMDRS{T}}
    version ∈ [:v1, :v2] || error("version must be v1 or v2")
    if version == :v1
        return _invfourier_mdrs_v1(kRvectors, Oᴿ, kpoints)
    else
        return _invfourier_mdrs_v2(kRvectors, Oᴿ, kpoints)
    end
end
