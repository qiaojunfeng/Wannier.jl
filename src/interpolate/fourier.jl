using LinearAlgebra

@doc raw"""
Fourier transform operator from k space to R space.
The so-called Wiger-Seitz (WS) interpolation.

```math
X_{mn}(\bf{R}) =
\frac{1}{N_{\bf{k}}} \sum_{\bf k} e^{-i{\bf k}{\mathbf{R}}} X_{mn}(\bf k)
```

O: operator matrix, n_wann * n_wann * n_kpts
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
Inverse Fourier transform operator from R space to k space.

```math
X_{mn}(\bf{k}) =
\sum_{\bf R} \frac{1}{N_{\bf R}} e^{i{\bf k}{\mathbf{R}}} X_{mn}(\bf R)
```

O_R: real space operator, i.e. in the frequency domain
kpoints: kpoints to be interpolated, 3 * n_kpts, can be nonuniform, fractional coordinates.
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
Fourier MDRS using plain loops.
Slower, but reproduce W90 behavior, and generate the same seedname_tb.dat file.

```math
X_{mn}(\bf R) =
\frac{1}{N_{\bf k}} \sum_{{\bf k}}e^{-i{\bf k}{\bf R}} X_{mn}(\bf k)
```
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
Faster version of Fourier MDRS, remove some for loops.
But need to expand R vectors to R+T, so different from seeedname_tb.dat.

```math
\widetilde{X}_{mn}({\tilde{\bf R}}) =
\sum_{{\bf R}} \frac{1}{N_{\bf R} {\cal N}_{mn{\bf R}}} X_{mn}({\bf R})
\sum_{j=1}^{{\cal N}_{mn{\bf R}}}   \delta_{{\tilde{\bf R}},{\bf R}+\mathbf{T}_{mn{\bf R}}^{(j)}}
```
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

```math
X_{mn}({\bf k}) =
\sum_{\bf R} \frac{1}{{\cal N}_{mn{\bf R}}} X_{mn}({\bf R})
\sum_{j=1}^{{\cal N}_{mn{\bf R}}} e^{i{\bf k}\cdot\left({\bf R}+\mathbf{T}_{mn{\bf R}}^{(j)}\right)}
```
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
```math
X_{mn}({\bf k}) =
\sum_{\tilde{\bf R}}e^{i{\bf k}{\tilde{\bf R}}} \widetilde{X}_{mn}({\tilde{\bf R}})
```
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
