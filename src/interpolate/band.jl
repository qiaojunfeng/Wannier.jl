import LinearAlgebra as LA
using FFTW
# import FINUFFT


"""
Get a matrix of kpoint coordinates from a kpath defined in win file.

kpath: fractional coordiantes.
n_points: number of kpoints in the first segment, remainning segments
    have the same density as 1st segment.
return kpoints in fractional coordinates.
"""
function get_kpath_points(kpath::Kpath{T}, n_points::Int, recip_lattice::AbstractMatrix{T}) where {T<:Real}
    # Use the kpath density of first segment to generate the following kpaths,
    # also need to take care of high symmetry kpoints at the start and end of each segment.
    n_seg = length(kpath)
    dk = 0.0
    ∑seg = 0.0

    # x axis value for plotting, cartesian length
    x = Vector{Float64}()

    # actual coordiantes, cartesian
    kpt = Matrix{Float64}(undef, 3, 0)

    # symmetry points
    n_symm = 0
    symm_idx = Vector{Int}()
    symm_label = Vector{String}()

    atol = 1e-7

    k1_prev = nothing
    k2_prev = nothing

    for i = 1:n_seg
        # a Pair: "L" => [0.5, 0.5, 0.5]
        (lab1, k1), (lab2, k2) = kpath[i]

        # to cartesian coordinates
        seg = recip_lattice * (k2 - k1)
        seg_norm = LA.norm(seg)

        if i == 1
            # kpath density
            dk = seg_norm / n_points
        end

        x_seg = collect(∑seg:dk:(∑seg + seg_norm))
        n_x_seg = length(x_seg)
        dk_seg = seg / seg_norm * dk

        # column vector * row vector = matrix
        kpt_seg = dk_seg * collect(1:n_x_seg)'
        kpt_seg .+= recip_lattice * k1

        if abs(∑seg + seg_norm - x_seg[end]) > atol
            # always include the last point
            push!(x_seg, ∑seg + seg_norm)
            k2_cart = recip_lattice * k2
            kpt_seg = hcat(kpt_seg, k2_cart)
        end

        if k2_prev !== nothing && all(isapprox.(k1, k2_prev; atol=atol))
            # remove repeated points
            popfirst!(x_seg)
            kpt_seg = kpt_seg[:, 2:end]
            n_symm += 1
        else
            n_symm += 2
            push!(symm_idx, length(x) + 1)
            push!(symm_label, lab1)
        end

        push!(symm_idx, length(x) + length(x_seg))
        push!(symm_label, lab2)

        append!(x, x_seg)
        kpt = hcat(kpt, kpt_seg)

        ∑seg += seg_norm

        k1_prev = k1
        k2_prev = k2
    end

    # to fractional
    kpt_frac = zeros(Float64, 3, length(x))
    kpt_frac = LA.inv(recip_lattice) * kpt

    return kpt_frac, x, symm_idx, symm_label
end

"""
nonuniform ifft for several kpoints.

    O_R: real space operator, i.e. in the frequency domain
    kpoints: kpoints in fractional coordinates, 3 x n_kpts
"""
function nuifft(O_R::Array{Complex{T}, 5}, kpoints::Matrix{T}) where {T<:Real}
    m, n, nx, ny, nz = size(O_R)

    nk = size(kpoints, 2)
    O_k = zeros(Complex{T}, m, n, nk)

    for ik = 1:nk
        for i = 1:nx, j = 1:ny, k = 1:nz
            kpt = kpoints[:, ik]
            fac = exp(im * 2π * LA.dot(kpt, [i-1;j-1;k-1]))
            O_k[:, :, ik] += fac * O_R[:, :, i, j, k]
        end
    end

    return O_k
end

"""
From kspace operator defined on a list of kpoints,
to a 5 dimensional array defined on x, y, z.
"""
function Ok_xyz(Ok::Array{T, 3}, xyz_k::Array{Int,3}) where {T<:Number}
    m, n, nk = size(Ok)
    nx, ny, nz = size(xyz_k)

    nk != nx * ny * nz && error("nk != nx * ny * nz")

    O_xyz = similar(Ok, m, n, nx, ny, nz)

    for i = 1:nx, j = 1:ny, k = 1:nz
        O_xyz[:, :, i, j, k] = Ok[:, :, xyz_k[i, j, k]]
    end

    return O_xyz
end


"""
interpolate band structure along a kpath
"""
function interpolate(model::Model{T}, kpoints::Matrix{T}) where {T<:Real}

    n_kx, n_ky, n_kz = model.kgrid
    k_xyz, xyz_k = get_kpoint_mappings(model.kpoints, model.kgrid)

    n_wann = model.n_wann

    H_k = zeros(Complex{T}, n_wann, n_wann, n_kx, n_ky, n_kz)
    # rotate
    for i = 1:n_kx, j = 1:n_ky, k = 1:n_kz
        ik = xyz_k[i, j, k]
        H_k[:, :, i, j, k] = model.A[:, :, ik]' * LA.Diagonal(model.E[:, ik]) * model.A[:, :, ik]
    end

    H_R = zeros(Complex{T}, n_wann, n_wann, n_kx, n_ky, n_kz)

    # bring to R space
    # for m = 1:data.n_wann
    #     for n = 1:data.n_wann
    #         ham_R[:,:,:,m,n] = FFTW.fft(ham_k[:,:,:,m,n], [1,2,3])
    #     end
    # end
    H_R .= FFTW.fft(H_k, [3, 4, 5])

    # default fftfreq(4, 1) = [0.0  0.25  -0.5  -0.25]
    # same as kmesh.pl, but if user use a different kgrid,
    # the results is wrong.
    model.kpoints[:, 1] ≉ zeros(T, 3) && error("kpoints[:, 0] ≉ zeros(3)")

    # tol = 1e-10
    # opts = FINUFFT.finufft_default_opts()
    # opts.modeord = 1
    # ham_kpath = zeros(ComplexF64, kpaths_num_kpts_tot, data.n_wann, data.n_wann)
    # for m = 1:data.n_wann, n = 1:data.n_wann
    #     # TODO: nufft3d2many
    #     ham_kpath[:, m, n] = FINUFFT.nufft3d2(
    #         kx .* 2 * pi, ky .* 2 * pi, kz .* 2 * pi,
    #         1, tol, ham_R[:, :, :, m, n], opts)
    #     # ham_kpath[:,m,n] = nuifft(
    #     #     kx, ky, kz, ham_R[:,:,:,m,n])
    # end
    # ham_kpath ./= data.num_kpts

    H_kpath = nuifft(H_R, kpoints)
    H_kpath ./= n_kx * n_ky * n_kz

    n_kpath_points = size(kpoints, 2)
    # diagonalize
    E_kpath = zeros(T, n_wann, n_kpath_points)
    for ik = 1:n_kpath_points
        H = H_kpath[:, :, ik]
        # @assert LA.ishermitian(H) H
        @warn LA.norm(H - H') ik
        H = 0.5 * (H + H')
        F = LA.eigen(H)
        E_kpath[:, ik] = real.(F.values)
    end

    return E_kpath
end
