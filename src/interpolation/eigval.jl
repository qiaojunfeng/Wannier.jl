import LinearAlgebra as LA
import FFTW
import FINUFFT
import Serialization


"""
kpts: in reduced coordinates, 3 * num_kpts
"""
function get_kpts_mapping(kpts::Matrix{Float64}, kpts_size::Vector{Int})
    num_kpts = prod(kpts_size)
    nk1, nk2, nk3 = kpts_size
    dk1, dk2, dk3 = 1 / nk1, 1 / nk2, 1 / nk3
    kpts_int = kpts ./ [dk1; dk2; dk3]
    kpts_int = round.(Int, kpts_int)
    for ik = 1:num_kpts
        kpts_int[1, ik] = mod(kpts_int[1, ik], 0:nk1-1) + 1
        kpts_int[2, ik] = mod(kpts_int[2, ik], 0:nk2-1) + 1
        kpts_int[3, ik] = mod(kpts_int[3, ik], 0:nk3-1) + 1
    end

    k2ijk = zeros(Int, num_kpts, 3)
    ijk2k = zeros(Int, nk1, nk2, nk3)
    for ik = 1:num_kpts
        k2ijk[ik, :] = kpts_int[:, ik]
        ijk2k[kpts_int[:, ik]...] = ik
    end

    return k2ijk, ijk2k
end

"""
kpath_startend: reduced coordiantes, 1st dim: vector x,y,z, 2nd dim: start,end, 3rd dim: num_segments
xyz_cart: whether return Cartesian coordiantes
"""
function generate_kpath_points(kpath_startend::Array{Float64,3}, kpath_label::Array{String,2}, bands_num_points::Int, recip_cell::Matrix{Float64}; xyz_cart::Bool=true)
    # Use the kpath density of first segment to generate the following kpaths,
    # also need to take care of high symmetry kpoints at the start and end of each segment.
    num_segments = size(kpath_startend, 3) # number of segments
    dk = 0
    sum_segments = 0.0
    kpaths = Vector{Float64}() # x axis value for plotting
    # actual coordiantes
    kx = Vector{Float64}()
    ky = Vector{Float64}()
    kz = Vector{Float64}()
    # symmetry points
    num_symm_points = 0
    tol = 1e-7
    symm_points = Vector{Int}()
    symm_points_label = Vector{String}()
    k1_prev = nothing
    k2_prev = nothing
    for i = 1:num_segments
        k1 = kpath_startend[:, 1, i]
        k2 = kpath_startend[:, 2, i]
        seg = recip_cell * (k2 - k1)
        seg_norm = LA.norm(seg)
        if i == 1
            # kpath density
            dk = seg_norm / bands_num_points
        end
        x = collect(sum_segments:dk:(sum_segments+seg_norm))
        dk_vec = seg / seg_norm * dk
        # column vector * row vector = matrix
        ik_vec = reshape(dk_vec, (3, 1)) .* reshape(1:length(x), (1, length(x)))
        ik_vec .+= recip_cell * k1
        ikx = ik_vec[1, :] # isa Vector, can use push! on it
        iky = ik_vec[2, :]
        ikz = ik_vec[3, :]
        if abs(sum_segments + seg_norm - x[end]) > tol
            # always include the last point
            push!(x, sum_segments + seg_norm)
            k2_cart = recip_cell * k2
            push!(ikx, k2_cart[1])
            push!(iky, k2_cart[2])
            push!(ikz, k2_cart[3])
        end
        if k2_prev !== nothing && all(isapprox.(k1, k2_prev; atol=tol))
            # remove repeated points
            popfirst!(x)
            popfirst!(ikx)
            popfirst!(iky)
            popfirst!(ikz)
            num_symm_points += 1
        else
            num_symm_points += 2
            push!(symm_points, length(kpaths) + 1)
            push!(symm_points_label, kpath_label[1, i])
        end
        push!(symm_points, length(kpaths) + length(x))
        push!(symm_points_label, kpath_label[2, i])
        append!(kpaths, x)
        append!(kx, ikx)
        append!(ky, iky)
        append!(kz, ikz)
        sum_segments += seg_norm
        k1_prev = k1
        k2_prev = k2
    end

    if !xyz_cart
        k = zeros(Float64, 3, length(kx))
        k[1, :] = kx
        k[2, :] = ky
        k[3, :] = kz
        k = LA.inv(recip_cell) * k
        kx = k[1, :]
        ky = k[2, :]
        kz = k[3, :]
    end
    return kpaths, kx, ky, kz, symm_points, symm_points_label
end

function nuifft(kx, ky, kz, ham_R)
    nx, ny, nz = size(ham_R)
    nh = length(kx)
    h = zeros(ComplexF64, nh)
    for ih = 1:nh
        for i = 1:nx, j = 1:ny, k = 1:nz
            h[ih] += exp(im * 2pi * (kx[ih] * (i - 1) + ky[ih] * (j - 1) + kz[ih] * (k - 1))) * ham_R[i, j, k]
        end
    end
    return h
end

"""
interpolate band structure along a kpath
"""
function interpolate(params::InputParams)

    data = read_seedname(params.seed_name; amn=true, mmn=false, eig=true)

    k2ijk, ijk2k = get_kpts_mapping(data.kpts, data.kpts_size)
    nk1, nk2, nk3 = data.kpts_size

    ham_k = zeros(ComplexF64, nk1, nk2, nk3, data.num_wann, data.num_wann)
    ham_R = zeros(ComplexF64, nk1, nk2, nk3, data.num_wann, data.num_wann)
    # rotate
    for i = 1:nk1, j = 1:nk2, k = 1:nk3
        ik = ijk2k[i, j, k]
        ham_k[i, j, k, :, :] = data.amn[:, :, ik]' * LA.Diagonal(data.eig[:, ik]) * data.amn[:, :, ik]
    end
    Serialization.serialize("ham_k.jls", ham_k)

    # bring to R space
    # for m = 1:data.num_wann
    #     for n = 1:data.num_wann
    #         ham_R[:,:,:,m,n] = FFTW.fft(ham_k[:,:,:,m,n], [1,2,3])
    #     end
    # end
    ham_R[:, :, :, :, :] = FFTW.fft(ham_k[:, :, :, :, :], [1, 2, 3])

    # nonuniform data: kpath
    kpaths, kx, ky, kz, symm_points, symm_points_label = generate_kpath_points(
        params.kpath, params.kpath_label, params.bands_num_points, data.recip_cell; xyz_cart=false)
    kpaths_num_kpts_tot = length(kpaths)
    num_symm_points = length(symm_points)

    tol = 1e-10
    opts = FINUFFT.finufft_default_opts()
    opts.modeord = 1
    ham_kpath = zeros(ComplexF64, kpaths_num_kpts_tot, data.num_wann, data.num_wann)
    for m = 1:data.num_wann, n = 1:data.num_wann
        # TODO: nufft3d2many
        ham_kpath[:, m, n] = FINUFFT.nufft3d2(
            kx .* 2 * pi, ky .* 2 * pi, kz .* 2 * pi,
            1, tol, ham_R[:, :, :, m, n], opts)
        # ham_kpath[:,m,n] = nuifft(
        #     kx, ky, kz, ham_R[:,:,:,m,n])
    end
    ham_kpath ./= data.num_kpts

    # diagonalize
    energies = zeros(Float64, kpaths_num_kpts_tot, data.num_wann)
    for i = 1:kpaths_num_kpts_tot
        # @assert LA.ishermitian(ham_kpath[i,:,:]) ham_kpath[i,:,:]
        @warn LA.norm(ham_kpath[i, :, :] - ham_kpath[i, :, :]') kx[i] ky[i] kz[i]
        ham_kpath[i, :, :] = 0.5 * (ham_kpath[i, :, :] + ham_kpath[i, :, :]')
        F = LA.eigen(ham_kpath[i, :, :])
        energies[i, :] = real.(F.values)
    end

    # output
    bands = Bands(
        num_kpts=kpaths_num_kpts_tot,
        num_bands=data.num_wann,
        kpaths=kpaths,
        kpaths_coord=vcat(kx', ky', kz'),
        energies=energies,
        num_symm_points=num_symm_points,
        symm_points=symm_points,
        symm_points_label=symm_points_label
    )

    write_wannier90_bands(bands, params.seed_name)
end
