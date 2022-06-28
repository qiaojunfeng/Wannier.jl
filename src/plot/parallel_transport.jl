import LinearAlgebra as LA
using CairoMakie


# Plot obstructions
function plot_surface_obstructions(p, suffix="")
    if p.N3 != 1
        phases = zeros(p.N1, p.N2, p.nwannier)
        for i = 1:p.N1, j = 1:p.N2
            Obs = normalize_matrix(overlap_A([i, j, p.N3], [i, j, 1], p)) #rotation at image point
            phases[i, j, :] = sort(imag(log.(eigvals(Obs))))
        end
        figure()
        xx = [p.t1[i] for i = 1:p.N1, j = 1:p.N2]
        yy = [p.t2[j] for i = 1:p.N1, j = 1:p.N2]
        for n = 1:p.nwannier
            plot_surface(xx, yy, phases[:, :, n], rstride=1, cstride=1)
        end
        savefig("wannierize$filename.$suffix.pdf")
        close()
    end
end


function plot_obstruction(model::Model{T}, A::Array{Complex{T},3}, obstruction::Obstruction{Complex{T}}, filename_prefix::String="obstruction") where {T<:Real}

    n_kx, n_ky, n_kz = model.kgrid
    k_xyz, xyz_k = get_kpoint_mappings(model.kpoints, model.kgrid)

    # Initialize arrays
    Ωright = zeros(n_kx, n_ky)
    Ωup = zeros(n_kx, n_ky)

    # Fill arrays
    imax = 1
    jmax = 1
    Ωmax = 0.0

    for i = 1:n_kx, j = 1:n_ky
        right = i % n_kx + 1
        up = j % n_ky + 1

        ik_right = xyz_k[right, j, 1]
        ik = xyz_k[i, j, 1]
        ik_up = xyz_k[i, up, 1]

        Mright = overlap(model.M, model.bvectors.kpb_k, ik_right, ik)
        Mup = overlap(model.M, model.bvectors.kpb_k, ik_up, ik)

        Ωright[i, j] = n_kx * LA.norm(A[:, :, ik_right] - Mright * A[:, :, ik])
        Ωup[i, j] = n_ky * LA.norm(A[:, :, ik_up] - Mup * A[:, :, ik])

        if Ωright[i, j] + Ωup[i, j] > Ωmax
            Ωmax = Ωright[i, j] + Ωup[i, j]
            imax = i
            jmax = j
        end
    end

    println("Ωmax = $Ωmax  [i,j] = $([imax,jmax])")

    fig = Figure()

    # display(Ωright)

    ax1 = Axis(fig[1, 1], aspect=DataAspect(), xlabel=L"$k_1$", ylabel=L"$k_2$",
     title=L"$\Omega_{right}$",
    )
    hm1 = heatmap!(ax1, Ωright; colormap=:binary)
    ax1.xticks = 1:n_kx
    ax1.yticks = 1:n_ky

    ax2 = Axis(fig[1, 2], aspect=DataAspect(), xlabel=L"$k_1$", ylabel=L"$k_2$", title=L"$\Omega_{up}$")
    hm2 = heatmap!(ax2, Ωup; colormap=:binary)
    ax2.xticks = 1:n_kx
    ax2.yticks = 1:n_ky

    colorrange = (min(minimum(Ωright), minimum(Ωup)), max(maximum(Ωright), maximum(Ωup)))
    Colorbar(fig[1, 3], colorrange=colorrange, colormap=:binary)

    # supertitle
    Label(fig[0, :], "Regularity of Bloch frame (finite difference with neighbor)")

    filename = "$filename_prefix.png"
    save(filename, fig; px_per_unit=3)
    println("Saved to $filename")

    return 

    plot_surface_obstructions(p, "_1_none")
    plot_surface_obstructions(p, "_2_corners")
    plot_surface_obstructions(p, "_3_edges")
    plot_surface_obstructions(p, "_4_surface")


    if n_kz == 1
        plot_3d = true #false
    else
        plot_3d = false
    end
    if plot_3d

        Uint = interp.Uint
        Obs_array = interp.Obs_array

        figure()
        for i in [div(n_ky, 3) div(n_ky, 3) * 2 n_ky] #5:5:n_ky #[33 66 100]
            plot3D(real(Uint[[1:n_ky; 1], i, 1, 1]), imag(Uint[[1:n_ky; 1], i, 1, 1]), real(Uint[[1:n_ky; 1], i, 2, 1]))#,"-k")
        end
        xlabel(L"$\mathrm{Re}(\mathbf{U}_{11})$")
        ylabel(L"$\mathrm{Im}(\mathbf{U}_{11})$")
        zlabel(L"$\mathrm{Re}(\mathbf{U}_{21})$")
        title(L"Interpolation of the first column of $s\mapsto\mathbf{U}(s)$ to $\mathbf{e}_1$")

        figure()
        for i in [div(n_ky, 3) div(n_ky, 3) * 2 n_ky] #5:5:n_ky #[33 66 100]
            plot3D(real(Uint[[1:n_ky; 1], i, 1, 1]), imag(Uint[[1:n_ky; 1], i, 1, 1]), imag(Uint[[1:n_ky; 1], i, 2, 1]), "-k")
        end
        xlabel(L"$\mathrm{Re}(\mathbf{U}_{11})$")
        ylabel(L"$\mathrm{Im}(\mathbf{U}_{11})$")
        zlabel(L"$\mathrm{Im}(\mathbf{U}_{21})$")
        title(L"Interpolation of the first column of $s\mapsto\mathbf{U}(s)$ to $\mathbf{e}_1$")

        eigs = zeros(ComplexF64, n_ky, 2)
        log_eig = zeros(n_ky, 2)
        c = zeros(2)

        for i in 2:n_ky
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
        plot(range(0, step=1, length=n_ky), log_eig[:, 1], label=L"$\log(\lambda_1)$")
        plot(range(0, step=1, length=n_ky), log_eig[:, 2], label=L"$\log(\lambda_2)$")
        xlabel(L"$k_1$")
        legend()
        title(L"Logarithm of the eigenvalues of $\mathbf{Obs}$")

        figure()
        plot3D(range(0, step=1, length=n_ky), cos.(log_eig[:, 1]), sin.(log_eig[:, 1]), label=L"$\lambda_1$")
        plot3D(range(0, step=1, length=n_ky), cos.(log_eig[:, 2]), sin.(log_eig[:, 2]), label=L"$\lambda_2$")
        legend()
        xlabel(L"$k_1$")
        ylabel(L"$\mathrm{Re}(\lambda_1)$")
        zlabel(L"$\mathrm{Im}(\lambda_2)$")
        title(L"Eigenvalues of $\mathbf{Obs}$")
    end

    plot_obs = false
    if (plot_obs)
        eig_obs = zeros(ComplexF64, n_kx, n_wann)
        eig_obs_i = zeros(ComplexF64, n_kx, n_wann)
        eig_obs_j = zeros(ComplexF64, n_kx, n_wann)

        for i = 1:n_kx
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


function plot_Bloch_frame_slice(p, A0, A0_init)
    if (log_interp & !p.map)
        suffix = "log"
    else
        suffix = "parallel_transport"
    end

    if n_kz > 2
        N3_list = [1, div(n_kz, 2), n_kz]
    else
        N3_list = [1]
    end
    init_frame = deepcopy(A0)
    for i = 1:n_kx, j = 1:n_ky, k = 1:n_kz
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
        savefig("Bloch_frame_$(p.filename)_sum_$(n_kx)_init.png")

        figure(figsize=(6.5, 5))
        matshow(real(init_frame[1, 1, :, :, slice] - init_frame[2, 2, :, :, slice]), false)
        colorbar()
        cur_axes = gca()
        cur_axes[:axis]("off")
        #title("Bloch frame slice before algorithm (A0), N3 = $slice")
        savefig("Bloch_frame_$(p.filename)_diff_$(n_kx)_init.png")

    end

    smooth_frame = deepcopy(A0)
    for i = 1:n_kx, j = 1:n_ky, k = 1:n_kz
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
        savefig("Bloch_frame_$(p.filename)_sum_$(n_kx)_$suffix.png")

        figure(figsize=(6.5, 5))
        matshow(real(smooth_frame[1, 1, :, :, slice] - smooth_frame[2, 2, :, :, slice]), false)
        colorbar()
        cur_axes = gca()
        cur_axes[:axis]("off")
        #title("Bloch frame slice after algorithm, N3 = $slice")
        savefig("Bloch_frame_$(p.filename)_diff_$(n_kx)_$suffix.png")
    end


end


#using Plots
#println("Creating animation of the obstruction")
#anim = @animate for i=1:2*p.N2-1
#               plot(real(Uint[:,abs(p.N2-i)+1,2,2]),imag(Uint[:,abs(p.N2-i)+1,2,2]),xlims=[-1,1],ylims=[-1,1])
#end

#gif(anim,"obs.gif",fps=25)
#close("all")
#Ucol1 = zeros(size(Uint)[1],size(Uint)[2],3)
#
#for i = 1:p.N1, j=1:p.N2
#	if Uint[i,j,2,1] != 0.0
#		phase = imag(log(im*Uint[i,j,2,1]/(abs(Uint[i,j,2,1]))))
#	else
#		phase = 0.0
#	end
#	array = Uint[i,j,:,1].*exp(-im*phase)
#
#	Ucol1[i,j,:] = [real(array[1]),imag(array[1]),imag(array[2])]
#end
