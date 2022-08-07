
function get_Hk(E::Matrix{T}, A::Array{U,3}) where {T<:Number,U<:Number}
    n_bands, n_wann, n_kpts = size(A)
    size(E) != (n_bands, n_kpts) && error("size(E) != (n_bands, n_kpts)")

    Hᵏ = zeros(U, n_wann, n_wann, n_kpts)
    for ik in 1:n_kpts
        Hᵏ[:, :, ik] = A[:, :, ik]' * Diagonal(E[:, ik]) * A[:, :, ik]
    end

    return Hᵏ
end

"""
interpolate band structure along a kpath
kpoints: interpolated kpoints in fractional coordinates, 3 x n_kpts, can be nonuniform.
"""
function interpolate(model::InterpolationModel{T}, kpoints::Matrix{T}) where {T<:Real}
    # n_wann x n_wann x n_kpts
    Hᵏ = get_Hk(model.E, model.A)

    # H_R = zeros(Complex{T}, n_wann, n_wann, n_kx, n_ky, n_kz)
    # bring to R space
    # for m = 1:n_wann
    #     for n = 1:n_wann
    #         H_R[m, n, :, :, :] = FFTW.fft(H_k[m, n, :, :, :], [3, 4, 5])
    #     end
    # end
    # H_R .= FFTW.fft(H_k, [3, 4, 5])

    # n_wann x n_wann x n_rvecs
    Hᴿ = fourier(model.kRvectors, Hᵏ)

    # default fftfreq(4, 1) = [0.0  0.25  -0.5  -0.25]
    # same as kmesh.pl, but if user use a different kgrid,
    # the results is wrong.
    # model.kpoints[:, 1] ≉ zeros(T, 3) && error("kpoints[:, 0] ≉ zeros(3)")

    Hᵏ_path = invfourier(model.kRvectors, Hᴿ, kpoints)

    # diagonalize
    Eᵏ_path = zeros(T, model.n_wann, size(kpoints, 2))
    for ik in axes(Eᵏ_path, 2)
        H = Hᵏ_path[:, :, ik]
        # @assert ishermitian(H) norm(H - H')
        @assert norm(H - H') < 1e-10
        # H = 0.5 * (H + H')
        ϵ, v = eigen(H)
        Eᵏ_path[:, ik] = real.(ϵ)
    end

    return Eᵏ_path
end

"""
Interpolate band structure along kpath.
"""
function interpolate(model::InterpolationModel, kpi::KPathInterpolant)
    # to Matrix
    kpoints = zeros(Float64, 3, length(kpi))
    for i in axes(kpoints, 2)
        kpoints[:, i] = kpi[i]
    end

    return interpolate(model, kpoints)
end

"""
Interpolate band structure along kpath.
"""
function interpolate(model::InterpolationModel)
    kpi = interpolate_w90(model.kpath, 100)
    return interpolate(model, kpi)
end
