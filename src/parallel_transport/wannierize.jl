include("interpolation.jl")
include("wannierize_utils.jl")

using PyPlot
using LinearAlgebra
#using Interpolations

## Assumptions: the kpoints are contained in a NxNxN cartesian grid, the neighbor list must contain the six cartesian neighbors

## N1xN2xN3 grid, filename.mmn must contain the overlaps
## nbeg and nend specify the window to wannierize
## Input Mmn file is nband x nband x nkpt x nntot
## Output is (nwannier = nend - nbeg + 1) x nband x nkpt, padded with zeros
function make_wannier(p, method)
    t1 = collect(0:p.N1-1) / p.N1
    t2 = collect(0:p.N2-1) / p.N2
    t3 = collect(0:p.N3-1) / p.N3
    Ntot = p.N1 * p.N2 * p.N3

    if lowercase(method) == "parallel transport"
        p.logMethod = false
    elseif lowercase(method) == "log interpolation"
        p.logMethod = true
    else
        println("method of interpolation '$method' not recognized")
        return
    end


    A0_init = deepcopy(p.A)
    M0 = deepcopy(p.M)
    A = p.A
    M = p.M
    neighbors = p.neighbors

    fill!(A, NaN) #protection: A must be filled by the algorithm
    println("Filling (k,0,0)")
    A[:, :, :, 1, 1] = propagate(Matrix(1.0I, p.nwannier, p.nwannier), [p.ijk_to_K[i, 1, 1] for i = 1:p.N1], p)


    # compute obstruction matrix
    Obs = normalize_matrix(overlap_A([p.N1, 1, 1], [1, 1, 1], p))
    println("Obstruction matrix = ")
    println(Obs)
    dV = eigen(Obs)
    d = dV.values
    V = dV.vectors
    logd = log.(d)
    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i = 1:p.nwannier
        if imag(logd[i]) < -pi + 0.01
            logd[i] = logd[i] + 2pi * im
        end
    end
    println("log(d) = $(logd)")

    # and pull it back
    for i = 1:p.N1
        A[:, :, i, 1, 1] = A[:, :, i, 1, 1] * V * diagm(0 => exp.(t1[i] * logd)) * V'
    end

    println("Filling (k1,k2,0)")
    for i = 1:p.N1
        A[:, :, i, :, 1] = propagate(A[:, :, i, 1, 1], [p.ijk_to_K[i, j, 1] for j = 1:p.N2], p)
    end
    A0 = deepcopy(A)

    # corner obstruction
    Obs = normalize_matrix(overlap_A([1, p.N2, 1], [1, 1, 1], p))
    dV = eigen(Obs)
    d = dV.values
    V = dV.vectors
    logd = log.(d)
    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i = 1:p.nwannier
        if imag(logd[i]) < -pi + 0.01
            logd[i] = logd[i] + 2pi * im
        end
    end
    # pull it back
    for i = 1:p.N1
        for j = 1:p.N2
            A[:, :, i, j, 1] = A[:, :, i, j, 1] * V * diagm(0 => exp.(t2[j] * logd)) * V'
        end
    end

    # Pull back the line obstruction
    phases = zeros(p.N1, p.nwannier)
    Obs_array = zeros(ComplexF64, p.N1, p.nwannier, p.nwannier)
    detObs = zeros(ComplexF64, p.N1)
    eigs = zeros(ComplexF64, p.N1, p.nwannier)
    for i = 1:p.N1
        Obs_array[i, :, :] = normalize_matrix(overlap_A([i, p.N2, 1], [i, 1, 1], p))
        detObs[i] = det(Obs_array[i, :, :])
    end

    # Find a continuous log of the determinant
    logDet = imag(log.(detObs))
    for i = 2:p.N1
        kmin = argmin([abs(logDet[i] + 2 * pi * k - logDet[i-1]) for k in -1:1])
        logDet[i] = logDet[i] + (kmin - 2) * 2 * pi
    end
    for i = 1:p.N1
        Obs_array[i, :, :] = exp(-im * logDet[i] / p.nwannier) * Obs_array[i, :, :]
        eigs[i, :] = eigvals(Obs_array[i, :, :])
    end

    # Interpolate the line obstruction
    Uint = zeros(ComplexF64, p.N1, p.N2, p.nwannier, p.nwannier)

    if !p.logMethod
        Uint = matrixTransport(Obs_array, t2)
    end

    for i = 1:p.N1
        Obs = Obs_array[i, :, :]
        dV = eigen(Obs)
        d = dV.values
        V = dV.vectors
        if p.logMethod #& (!p.map)
            for j = 1:p.N2
                Uint[i, j, :, :] = powm(Obs, t2[j])
            end
        end

        for j = 1:p.N2
            A[:, :, i, j, 1] = A[:, :, i, j, 1] * exp(im * logDet[i] * t2[j] / 2)
            A[:, :, i, j, 1] = A[:, :, i, j, 1] * Uint[i, j, :, :]
        end
    end


    # Propagate along the third dimension
    println("Filling (k1,k2,k3)")
    for i = 1:p.N1, j = 1:p.N2
        A[:, :, i, j, :] = propagate(A[:, :, i, j, 1], [p.ijk_to_K[i, j, k] for k = 1:p.N3], p)
    end

    # Fix corner
    Obs = normalize_matrix(overlap_A([1, 1, p.N3], [1, 1, 1], p))
    dV = eigen(Obs)
    d = dV.values
    V = dV.vectors
    logd = imag(log.(d))
    #println("Obstruction matrix  - Id = $(norm(Obs-I)))")

    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i = 1:p.nwannier
        if imag(logd[i]) < -pi + 0.01
            logd[i] = logd[i] + 2pi * im
        end
    end

    for k = 1:p.N3
        # fixer = powm(Obs,t3[k])
        fixer = V * diagm(0 => exp.(t3[k] * logd)) * V'
        for i = 1:p.N1, j = 1:p.N2
            A[:, :, i, j, k] = A[:, :, i, j, k] * fixer
        end
    end

    # Fix first edge
    Obs_array_i = zeros(ComplexF64, p.N1, p.nwannier, p.nwannier)
    Uint_ik = zeros(ComplexF64, p.N1, p.N1, p.nwannier, p.nwannier)
    for i = 1:p.N1
        Obs_array_i[i, :, :] = normalize_matrix(overlap_A([i, 1, p.N3], [i, 1, 1], p))
    end
    if !p.logMethod
        Uint_ik = matrixTransport(Obs_array_i, t3)
    end
    for i = 1:p.N1
        for k = 1:p.N3
            if p.logMethod
                fixer = powm(Obs_array_i[i, :, :], t3[k])
                for j = 1:p.N3
                    A[:, :, i, j, k] = A[:, :, i, j, k] * fixer
                end
            else
                for j = 1:p.N3
                    A[:, :, i, j, k] = A[:, :, i, j, k] * Uint_ik[i, k, :, :]
                end
            end
        end
    end

    # Fix second edge
    Obs_array_j = zeros(ComplexF64, p.N2, p.nwannier, p.nwannier)
    Uint_jk = zeros(ComplexF64, p.N1, p.N1, p.nwannier, p.nwannier)
    for j = 1:p.N2
        Obs_array_j[j, :, :] = normalize_matrix(overlap_A([1, j, p.N3], [1, j, 1], p))
    end
    if !p.logMethod
        Uint_jk = matrixTransport(Obs_array_j, t3)
    end
    for j = 1:p.N2
        for k = 1:p.N3
            if p.logMethod
                fixer = powm(Obs_array_j[j, :, :], t3[k])
                for i = 1:p.N1
                    A[:, :, i, j, k] = A[:, :, i, j, k] * fixer
                end
            elseif !p.logMethod
                for i = 1:p.N1
                    A[:, :, i, j, k] = A[:, :, i, j, k] * Uint_jk[j, k, :, :]
                end
            end
        end
    end

    # Fix whole surface
    for i = 1:p.N1, j = 1:p.N2
        Obs = normalize_matrix(overlap_A([i, j, p.N3], [i, j, 1], p))
        for k = 1:p.N3
            A[:, :, i, j, k] = A[:, :, i, j, k] * powm(Obs, t3[k])
        end
    end


    p.A = A
    interp = InterpResults(Obs_array, Obs_array_i, Obs_array_j, Uint, Uint_ik)

    return (p, interp)
end

function plot_results(p, interp)
    #Interpolation method
    if (p.logMethod)
        suffix = "log"
    else
        suffix = "parallel_transport"
    end

    #Initialize arrays
    omegaright = zeros(p.N1, p.N2)
    omegaup = zeros(p.N1, p.N2)
    Aright1 = zeros(ComplexF64, p.N1, p.N2)
    Aright2 = zeros(ComplexF64, p.N1, p.N2)
    detA = zeros(ComplexF64, p.N1, p.N2)

    #Fill arrays
    imax = 1
    jmax = 1
    omegamax = 0.0
    for i = 1:p.N1, j = 1:p.N2
        right = i == p.N1 ? 1 : i + 1
        up = j == p.N2 ? 1 : j + 1
        Krightj = p.ijk_to_K[right, j, 1]
        Kij = p.ijk_to_K[i, j, 1]
        Kiup = p.ijk_to_K[i, up, 1]
        omegaright[i, j] = (p.N1) * norm(p.A[:, :, right, j, 1] - overlap(Krightj, Kij, p) * p.A[:, :, i, j, 1])
        omegaup[i, j] = (p.N2) * norm(p.A[:, :, i, up, 1] - overlap(Kiup, Kij, p) * p.A[:, :, i, j, 1])
        if (omegaright[i, j] + omegaup[i, j] > omegamax)
            omegamax = omegaright[i, j] + omegaup[i, j]
            imax = i
            jmax = j
        end

    end
    println("omegamax = $omegamax, [i,j] = $([imax,jmax])")
    if p.N1 > 10
        nb = 10
    else
        nb = 1
    end
    fig_size = (4, 3)

    figure(figsize=fig_size)
    ax = PyPlot.axes()
    matshow(omegaright, false, cmap=ColorMap("binary"), origin="lower")
    xticks(0:(div(p.N1, nb)):p.N1, 0:(1.0/nb):1)
    yticks(0:(div(p.N1, nb)):p.N1, 0:(1.0/nb):1)
    ax[:xaxis][:set_ticks_position]("bottom")
    colorbar()
    xlabel(L"$k_1$")
    ylabel(L"$k_2$")
    #cur_axes = gca()
    #cur_axes[:axis]("off")
    savefig("omega_right_$(p.filename)_$suffix.png")
    title("Regularity of Bloch frame (finite difference with right neighbor)")

    figure(figsize=fig_size)
    ax = PyPlot.axes()
    matshow(omegaup, false, cmap=ColorMap("binary"), origin="lower")
    xticks(0:(div(p.N1, nb)):p.N1, 0:(1.0/nb):1)
    yticks(0:(div(p.N1, nb)):p.N1, 0:(1.0/nb):1)
    ax[:xaxis][:set_ticks_position]("bottom")
    colorbar()
    xlabel(L"$k_1$")
    ylabel(L"$k_2$")
    #cur_axes = gca()
    #cur_axes[:axis]("off")
    savefig("omega_up_$(p.filename)_$suffix.png")
    title("Regularity of Bloch frame (finite difference with up neighbor)")

    plot_surface_obstructions(p, "_1_none")
    plot_surface_obstructions(p, "_2_corners")
    plot_surface_obstructions(p, "_3_edges")
    plot_surface_obstructions(p, "_4_surface")


    if p.N3 == 1
        plot_3d = true #false
    else
        plot_3d = false
    end
    if plot_3d

        Uint = interp.Uint
        Obs_array = interp.Obs_array

        figure()
        for i in [div(p.N2, 3) div(p.N2, 3) * 2 p.N2] #5:5:p.N2 #[33 66 100]
            plot3D(real(Uint[[1:p.N2; 1], i, 1, 1]), imag(Uint[[1:p.N2; 1], i, 1, 1]), real(Uint[[1:p.N2; 1], i, 2, 1]))#,"-k")
        end
        xlabel(L"$\mathrm{Re}(\mathbf{U}_{11})$")
        ylabel(L"$\mathrm{Im}(\mathbf{U}_{11})$")
        zlabel(L"$\mathrm{Re}(\mathbf{U}_{21})$")
        title(L"Interpolation of the first column of $s\mapsto\mathbf{U}(s)$ to $\mathbf{e}_1$")

        figure()
        for i in [div(p.N2, 3) div(p.N2, 3) * 2 p.N2] #5:5:p.N2 #[33 66 100]
            plot3D(real(Uint[[1:p.N2; 1], i, 1, 1]), imag(Uint[[1:p.N2; 1], i, 1, 1]), imag(Uint[[1:p.N2; 1], i, 2, 1]), "-k")
        end
        xlabel(L"$\mathrm{Re}(\mathbf{U}_{11})$")
        ylabel(L"$\mathrm{Im}(\mathbf{U}_{11})$")
        zlabel(L"$\mathrm{Im}(\mathbf{U}_{21})$")
        title(L"Interpolation of the first column of $s\mapsto\mathbf{U}(s)$ to $\mathbf{e}_1$")

        eigs = zeros(ComplexF64, p.N2, 2)
        log_eig = zeros(p.N2, 2)
        c = zeros(2)

        for i in 2:p.N2
            eigs[i, :], w = eigen(Obs_array[i, :, :])
            log_eig[i, :] = imag(log.(eigs[i, :]))
            log_eig[i, :] += c[:] * 2π
            if log_eig[i-1, 1] - log_eig[i, 1] > 1.8π
                c[1] += 1
                if log_eig[i-1, 1] - log_eig[i, 1] > 3.8π
                    c[1] += 1
                end
            else
                c[1] = 0
            end
            if log_eig[i, 2] - log_eig[i-1, 2] > 1.8π
                c[2] += -1
                if log_eig[i, 2] - log_eig[i-1, 2] > 3.8π
                    c[2] += -1
                end
            else
                c[2] = 0
            end

            log_eig[i, :] += c[:] * 2π
        end
        figure()
        plot(range(0, step=1, length=p.N2), log_eig[:, 1], label=L"$\log(\lambda_1)$")
        plot(range(0, step=1, length=p.N2), log_eig[:, 2], label=L"$\log(\lambda_2)$")
        xlabel(L"$k_1$")
        legend()
        title(L"Logarithm of the eigenvalues of $\mathbf{Obs}$")

        figure()
        plot3D(range(0, step=1, length=p.N2), cos.(log_eig[:, 1]), sin.(log_eig[:, 1]), label=L"$\lambda_1$")
        plot3D(range(0, step=1, length=p.N2), cos.(log_eig[:, 2]), sin.(log_eig[:, 2]), label=L"$\lambda_2$")
        legend()
        xlabel(L"$k_1$")
        ylabel(L"$\mathrm{Re}(\lambda_1)$")
        zlabel(L"$\mathrm{Im}(\lambda_2)$")
        title(L"Eigenvalues of $\mathbf{Obs}$")
    end

    plot_obs = false
    if (plot_obs)
        eig_obs = zeros(ComplexF64, p.N1, p.nwannier)
        eig_obs_i = zeros(ComplexF64, p.N1, p.nwannier)
        eig_obs_j = zeros(ComplexF64, p.N1, p.nwannier)

        for i = 1:p.N1
            eig_obs[i, :] = eigen(Obs_array[i, :, :])[1]
            eig_obs_i[i, :] = eigen(Obs_array_i[i, :, :])[1]
            eig_obs_j[i, :] = eigen(Obs_array_j[i, :, :])[1]
        end

        figure()
        plot(real(eig_obs))
        figure()
        plot(real(eig_obs_i))
        figure()
        plot(real(eig_obs_j))
    end
end

function print_error(p)
    err_before = 0.0
    err_after = 0.0
    for i = 1:p.N1, j = 1:p.N2, k = 1:p.N3
        err_before += norm(normalize_matrix(overlap(p.ijk_to_K[i, j, k], p.ijk_to_K[(i%p.N1+1), j, k], p)) - I)^2
        err_after += norm(normalize_matrix(overlap_A([i, j, k], [(i % p.N1 + 1), j, k], p)) - I)^2
        err_before += norm(normalize_matrix(overlap(p.ijk_to_K[i, j, k], p.ijk_to_K[i, (j%p.N2)+1, k], p)) - I)^2
        err_after += norm(normalize_matrix(overlap_A([i, j, k], [i, (j % p.N2) + 1, k], p)) - I)^2
    end
    err_before = p.N1 * sqrt(1 / (p.N1 * p.N2 * p.N3) * err_before)
    err_after = p.N1 * sqrt(1 / (p.N1 * p.N2 * p.N3) * err_after)
    println("err before = $(err_before)")
    println("err after  = $(err_after)")
end

function plot_Bloch_frame_slice(p, A0, A0_init)
    if (p.logMethod & !p.map)
        suffix = "log"
    else
        suffix = "parallel_transport"
    end

    if p.N3 > 2
        N3_list = [1, div(p.N3, 2), p.N3]
    else
        N3_list = [1]
    end
    init_frame = deepcopy(A0)
    for i = 1:p.N1, j = 1:p.N2, k = 1:p.N3
        init_frame[:, :, i, j, k] = A0_init[:, :, i, j, k]' * A0[:, :, i, j, k] #*p.A[:,:,i,j,k]
        #init_frame[:,:,i,j,k] = normalize_matrix(init_frame[:,:,i,j,k])
    end

    for slice in N3_list
        figure(figsize=(6.5, 5))
        matshow(real(init_frame[1, 1, :, :, slice] + init_frame[2, 2, :, :, slice]), false)
        colorbar()
        cur_axes = gca()
        cur_axes[:axis]("off")
        #title("Bloch frame slice before algorithm (A0), N3 = $slice")
        savefig("Bloch_frame_$(p.filename)_sum_$(p.N1)_init.png")

        figure(figsize=(6.5, 5))
        matshow(real(init_frame[1, 1, :, :, slice] - init_frame[2, 2, :, :, slice]), false)
        colorbar()
        cur_axes = gca()
        cur_axes[:axis]("off")
        #title("Bloch frame slice before algorithm (A0), N3 = $slice")
        savefig("Bloch_frame_$(p.filename)_diff_$(p.N1)_init.png")

    end

    smooth_frame = deepcopy(A0)
    for i = 1:p.N1, j = 1:p.N2, k = 1:p.N3
        smooth_frame[:, :, i, j, k] = p.A[:, :, i, j, k]' * A0[:, :, i, j, k] #*p.A[:,:,i,j,k]
        smooth_frame[:, :, i, j, k] = normalize_matrix(smooth_frame[:, :, i, j, k])
    end

    for slice in N3_list
        figure(figsize=(6.5, 5))
        matshow(real(smooth_frame[1, 1, :, :, slice] + smooth_frame[2, 2, :, :, slice]), false)
        colorbar()
        cur_axes = gca()
        cur_axes[:axis]("off")
        #title("Bloch frame slice after algorithm, N3 = $slice")
        savefig("Bloch_frame_$(p.filename)_sum_$(p.N1)_$suffix.png")

        figure(figsize=(6.5, 5))
        matshow(real(smooth_frame[1, 1, :, :, slice] - smooth_frame[2, 2, :, :, slice]), false)
        colorbar()
        cur_axes = gca()
        cur_axes[:axis]("off")
        #title("Bloch frame slice after algorithm, N3 = $slice")
        savefig("Bloch_frame_$(p.filename)_diff_$(p.N1)_$suffix.png")
    end


end
