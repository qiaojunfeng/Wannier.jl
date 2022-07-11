using PyPlot
using FFTW

linspace(a, b, n) = range(a; stop=b, length=n)

close("all")

## Uncomment for production
PyPlot.rc("font"; family="serif")
PyPlot.rc("xtick"; labelsize="x-small")
PyPlot.rc("ytick"; labelsize="x-small")
PyPlot.rc("figure"; figsize=(4, 3))
PyPlot.rc("text"; usetex=false)

nband = p.nband
Ks = zeros(Int, 3, nband)
Ntot = p.N1 * p.N2 * p.N3
Ns = [p.N1, p.N2, p.N3]

let ind = 1
    for l1 in (-L1):L1
        for l2 in (-L2):L2
            for l3 in (-L3):L3
                Ks[:, ind] = [l1, l2, l3]
                ind += 1
            end
        end
    end
end

Wfour = zeros(ComplexF64, nwannier, (2L1 + 1) * N1, (2L2 + 1) * N2, (2L3 + 1) * N3)
H = zeros(ComplexF64, nwannier, nwannier, N1, N2, N3)

ijk_to_kind(i, j, k) = k + j * N3 + i * N2 * N3 + 1
for i in 0:(N1 - 1)
    for j in 0:(N2 - 1)
        for k in 0:(N3 - 1)
            ijk = [i; j; k]
            kpt = ijk ./ Ns
            kind = ijk_to_kind(i, j, k)
            eig_k = λ(kpt, Ks)
            perm = (sortperm(eig_k))
            H[:, :, i + 1, j + 1, k + 1] =
                A[:, :, i + 1, j + 1, k + 1]' *
                diagm(0 => sort(eig_k)) *
                A[:, :, i + 1, j + 1, k + 1]
            @assert A[:, :, i + 1, j + 1, k + 1]'A[:, :, i + 1, j + 1, k + 1] ≈ I
            for iwann in 1:nwannier
                for iband in 1:nband
                    Wfour[iwann, i + 1 + (Ks[1, perm[iband]] + L1) * N1, j + 1 + (Ks[2, perm[iband]] + L2) * N2, k + 1 + (Ks[3, perm[iband]] + L3) * N3] = A[
                        iband, iwann, i + 1, j + 1, k + 1
                    ]
                    # mark visually the max of loc
                    # if j==jmax && k==kmax
                    #     Wfour[iwann,
                    #           i+1+(Ks[1,perm[iband]]+L1)*N1,
                    #           j+1+(Ks[2,perm[iband]]+L2)*N2,
                    #           k+1+(Ks[3,perm[iband]]+L3)*N3] = 0
                    # end
                end
            end
        end
    end
end

Wreal = zeros(ComplexF64, nwannier, (2L1 + 1) * N1, (2L2 + 1) * N2, (2L3 + 1) * N3)
# TODO this is a huge hack because I'm too lazy to do it properly. It makes WF be real when they should be
if do_shift == false && corner == [0.0, 0.0, 0.0] && (N1 == N2 == 1)
    new_Wfour = circshift(Wfour, (0, 0, 0, N3 / 2))
    for iwann in 1:nwannier
        Wreal[iwann, :, :, :] = ifftshift(ifft(fftshift(new_Wfour[iwann, :, :, :])))
    end
else
    for iwann in 1:nwannier
        Wreal[iwann, :, :, :] = ifftshift(ifft(fftshift(Wfour[iwann, :, :, :])))
    end
end

# See https://julialang.org/blog/2016/02/iteration for the trick used here for a completely generic (and fast!) multidimensional algorithm
function four_interp(arr, fac=1)
    R = CartesianIndices(size(arr))
    Ns = [size(arr)...]
    arr_fine = zeros(ComplexF64, (2fac - 1) * Ns...)
    shift = CartesianIndex((fac - 1) * Ns...)
    # function barrier for type inference
    function fillshift!(dest, src, R, shift)
        for I in R
            dest[shift + I] = src[I]
        end
    end
    fillshift!(arr_fine, fftshift(fft(arr)), R, shift)
    return (2fac - 1)^length(size(arr)) * ifft(ifftshift(arr_fine))
end

function four_interp_ham(H, fac=1)
    return mapslices(arr -> four_interp(arr, fac), H; dims=3:length(size(H)))
end

# 1D
if dim == 1
    kspace = (0:(N3 - 1)) / N3 .+ k_red_to_real([0.0, 0.0, 0.0])[3]
    ξspace = linspace(-L3 + kspace[1], L3 + kspace[end], (2L3 + 1) * N3) #TODO not actually linspace, but not that important
    rspace = linspace(-pi * N3, pi * N3, (2L3 + 1) * N3)
    # plot the gauge
    figure()
    plot(kspace, sort(p.eigs; dims=2))
    ylim(0, 2)
    xlim(kspace[1], kspace[end])
    xlabel("k")
    ylabel("ɛk")
    tight_layout()
    savefig("bandplot.pdf")

    figure()
    markers = ("s", "x", "o", "d")
    markers = ("-x", "--s", "o", "d")
    for n in 2:2
        # for n=1:1
        gca()[:set_color_cycle](nothing)
        # for m=1:p.nband
        for m in 1:3
            plot(
                kspace,
                real.(A[m, n, 1, 1, :]),
                markers[n];
                label="A$m$n",
                markersize=4,
                markeredgewidth=2,
                linewidth=2,
            )
            plot(
                kspace,
                imag.(A[m, n, 1, 1, :]),
                markers[n];
                label="A$m$n",
                markersize=4,
                markeredgewidth=2,
                linewidth=2,
            )
            # plot(imag.(A[m,n,1,1,:]),markers[n], label="$m $n",markersize=8, markeredgewidth=2,linewidth=2)
        end
    end
    xlim(kspace[1], kspace[end])
    legend()
    xlabel("k")
    ylabel("Re(Amn)")
    tight_layout()
    savefig("gauge.pdf")

    N = N3
    Wfour = dropdims(Wfour; dims=(2, 3))
    Wreal = dropdims(Wreal; dims=(2, 3))
    H = dropdims(H; dims=(3, 4))

    # normalize here rather than try to figure out correct normalizations
    for n in 1:(p.nwannier)
        Wfour[n, :] /= sqrt(sum(abs2, Wfour[n, :]) * (ξspace[2] - ξspace[1]))
        Wreal[n, :] /= sqrt(sum(abs2, Wreal[n, :]) * (rspace[2] - rspace[1]))
    end

    figure()
    plot(ξspace, real(Wfour)[1, :], "-")
    # plot(ξspace,real(Wfour)[2,:],"--")
    plot(ξspace, imag(Wfour)[1, :], "--")
    # plot(ξspace,imag(Wfour)[2,:],"-.")
    xlim((-3, 3))
    legend(vcat(["Re(w$i)" for i in 1:1], ["Im(w$i)" for i in 1:1]))
    xlabel("ξ")
    ylabel("w(ξ)")
    tight_layout()
    savefig("wfour.pdf")
    # plot(ξspace,imag(Wfour))

    figure()
    semilogy(rspace, abs.(Wreal)[1, :])
    # semilogy(rspace,1./rspace.^2,"k")
    # legend(["w1","1/r²"])
    xlim(rspace[1] / 2, rspace[end] / 2)
    ylim(1e-10, 1)
    xlabel("r")
    ylabel("|w(r)|")
    tight_layout()
    savefig("wreal_semilog.pdf")
    # plot(rspace,imag(Wreal))

    figure()
    print(size(Wreal))
    plot(rspace, real(Wreal)[1, :], "-")
    plot(rspace, real(Wreal)[2, :], "--")
    # plot(rspace,real(Wreal)[3,:],"--")
    # plot(rspace,imag(Wreal)')
    legend(["w1", "w2"])
    xlabel("r")
    ylabel("w(r)")
    tight_layout()
    savefig("wreal.pdf")
    # plot(rspace,imag(Wreal))

    figure()
    fac = 10
    Nint = (2fac - 1) * N
    x = collect((0:(N - 1)) / N)
    xint = collect((0:(Nint - 1)) / Nint)

    if do_shift
        # This is a bit awkward because Fourier interpolation works on a (0:N-1)/N grid
        x += 1 / 2 / N
        xint += 1 / 2 / N # this is really N and not Nint

        ## Here's a test showing why Nint: Fourier interpolation on a staggered grid (at least as defined here) actually extrapolates a bit also
        # N = 10
        # x = ((0:N-1)/N+1/2/N)*2pi
        # f = sin.(x)
        # plot(x,f)
        # Nint = 5N
        # xint = ((0:Nint-1)/Nint+1/2/Nint)*2pi
        # xint = ((0:Nint-1)/Nint+1/2/N)*2pi
        # plot(xint,four_interp(f,3))

    end

    x .+= corner[3]
    xint .+= corner[3]

    Hint = four_interp_ham(H, fac)

    # not in reduced BZ
    λ_real(k, Ks) = dropdims(sum(abs2, (Ks .+ k); dims=1); dims=1)

    plot(xint, [sort(λ_real([0, 0, xint[i]], Ks))[1] for i in 1:Nint], "-k"; label="Exact")
    plot(xint, [sort(λ_real([0, 0, xint[i]], Ks))[2] for i in 1:Nint], "-k")
    # plot(xint,[sort(λ_real([0,0,xint[i]],Ks))[3] for i =1:Nint], "-k")
    # ylim(-.1,1.1)
    # plot(xint,[sort(λ([0,0,xint[i]],Ks))[3] for i =1:Nint], "-k")
    plot(
        x,
        [sort(real(eigen(H[:, :, i]).values))[m] for i in 1:N, m in 1],
        "x";
        label="Samples",
    )
    plot(x, [sort(real(eigen(H[:, :, i]).values))[m] for i in 1:N, m in 2], "x")
    # plot(x,[sort(real(eigen(H[:,:,i]).values))[m] for i=1:N, m=3], "o")
    gca()[:set_color_cycle](nothing)
    plot(
        xint,
        [sort(real(eigen(Hint[:, :, i]).values))[m] for i in 1:Nint, m in 1],
        "-";
        label="Interpolation",
    )
    plot(xint, [sort(real(eigen(Hint[:, :, i]).values))[m] for i in 1:Nint, m in 2], "-")
    # plot(xint,[sort(real(eigen(Hint[:,:,i]).values))[m] for i=1:Nint, m=3], "-")
    legend()
    xlabel("k")
    ylabel("ɛ")
    xlim([corner[3], corner[3] + 1]) #clip to actual BZ, not more. TODO wrap stuff above corner[3]+1 around
    tight_layout()
    savefig("interp.pdf")

    maxerr = maximum(
        abs.(
            [sort(real(eigen(Hint[:, :, i]).values))[m] for i in 1:Nint, m in 1] - [sort(λ_real([0, 0, xint[i]], Ks))[1] for i in 1:Nint]
        ),
    )
    println("max err $maxerr")
end

if dim == 2
    npt = 100
    kx = linspace(0, 1, npt) .+ k_red_to_real([0.0, 0.0, 0.0])[2]
    ky = linspace(0, 1, npt) .+ k_red_to_real([0.0, 0.0, 0.0])[3]
    eigs = zeros(npt, npt, nband)
    for i in 1:npt
        for j in 1:npt
            eigs[i, j, :] = sort(λ([0, kx[i], ky[j]], Ks))
        end
    end
    # plot_surface(kspace*kspace',kspace*kspace',eigs[:,:,1])

    if (false)
        for n in 1:min(p.nband - 1, 6)
            # plot_surface(kspace,kspace,eigs[:,:,n],cmap=ColorMap("coolwarm"),linewidth=0,rstride=1,cstride=1,antialiased=false)
            # plot_surface(kspace,kspace,eigs[:,:,n+1] - eigs[:,:,n],cmap=ColorMap("coolwarm"),linewidth=0,rstride=1,cstride=1,antialiased=false)
            figure()
            pcolormesh(kx, ky, log10.(eigs[:, :, n + 1] - eigs[:, :, n] + 1e-3))
            title("$n -> $(n+1)")
            # pcolormesh(log10.(eigs[:,:,n]))
            colorbar()
        end
    end

    kx = linspace(0, 1, N2) .+ k_red_to_real([0.0, 0.0, 0.0])[2]
    ky = linspace(0, 1, N3) .+ k_red_to_real([0.0, 0.0, 0.0])[3]
    figure()
    pcolormesh(kx, ky, loc[1, :, :]; cmap=ColorMap("gray_r"))
    xlabel("kx")
    ylabel("ky")
    colorbar()
    savefig("normgrad_$N3.pdf")
    @assert N2 == N3
    N = N2

    Wfour = dropdims(Wfour; dims=(2))
    Wreal = dropdims(Wreal; dims=(2))
    H = dropdims(H; dims=(3))

    ξ = (0:(N3 - 1)) / N3 .+ k_red_to_real([0.0, 0.0, 0.0])[3]
    ξspacex = linspace(-L2 - kx[1], L2 + kx[end], (2L2 + 1) * N3) #TODO not actually linspace, but not that important
    ξspacey = linspace(-L3 - ky[1], L3 + ky[end], (2L3 + 1) * N3) #TODO not actually linspace, but not that important

    for iwann in 1:(p.nwannier)
        # plot the WF in Fourier space
        # figure()
        # pcolormesh(ξspacex,ξspacey,real(Wfour[iwann,:,:]))
        # figure()
        # pcolormesh(ξspacex,ξspacey,imag(Wfour[iwann,:,:]))
        figure()
        pl = pcolormesh(
            ξspacex,
            ξspacey,
            abs.(Wfour[iwann, :, :]);
            cmap=ColorMap("gray_r"),
            linewidth=0,
            rasterized=true,
        ) #force rasterization to avoid "kilt" effect when saving
        xlim(-2, 2)
        ylim(-2, 2)
        xlabel(L"\xi_x")
        ylabel(L"\xi_y")
        colorbar()
        tight_layout()
        savefig("2D_abs_$iwann.pdf"; dpi=450)

        # # plot the gradient in Fourier space
        # arr = real(Wfour[iwann,:,:])
        # grad = zeros(size(arr))
        # for i=1:size(arr,1)-1
        #     for j=1:size(arr,2)-1
        #         grad[i,j] = sqrt((arr[i+1,j]-arr[i,j])^2 + (arr[i,j+1]-arr[i,j])^2)
        #     end
        # end
        # grad /= (ξspacex[2]-ξspacex[1])
        # figure()
        # pcolormesh(ξspacex,ξspacey, grad,cmap=ColorMap("gray_r"),linewidth=0,rasterized=true)
        # xlim(-2,2)
        # ylim(-2,2)
        # xlabel(L"\xi_x")
        # ylabel(L"\xi_y")
        # colorbar()
        # tight_layout()
        # savefig("2D_grad_$iwann.pdf",dpi=450)
        # println("Max gradient $iwann: $(maximum(grad))")

        # plot a slice
        figure()
        data = abs.(Wreal[iwann, div(end, 2), :])
        rspace = linspace(-pi * N3, pi * N3, (2L3 + 1) * N3)
        semilogy(rspace, data)
        xlim(rspace[1] / 2, rspace[end] / 2)
        semilogy(rspace, 1 ./ rspace .^ 2 / 100, "k")
        legend(["w$iwann", "1/r²"])
        xlabel("r")
        ylabel("|w(r)|")
        tight_layout()
        savefig("wreal_slice_$iwann.pdf")

        figure()
        rspace = linspace(-pi * N3, pi * N3, (2L3 + 1) * N3)
        beg = 160
        finish = 200
        pcolormesh(rspace, rspace, real(Wreal[iwann, :, :]); linewidth=0, rasterized=true)
        xlim(-10, 10)
        ylim(-10, 10)
        xlabel("x")
        ylabel("y")
        colorbar()
        tight_layout()
        savefig("wreal_2D_$iwann.pdf"; dpi=450)
    end

    fac = 2
    Nint = (2fac - 1) * N
    x = (0:(N - 1)) / N
    xint = (0:(Nint - 1)) / Nint

    Hint = four_interp_ham(H, fac)

    @assert do_shift == false #not implemented yet

    figure()

    #exact
    plot(
        xint,
        [sort(λ([0, xint[i], 0], Ks))[m] for i in 1:Nint, m in 1:1],
        "-k";
        label="Exact",
    )
    plot(xint, [sort(λ([0, xint[i], 0], Ks))[m] for i in 1:Nint, m in 1:nwannier], "-k")
    plot(
        xint .+ 1, [sort(λ([0, 0, xint[i]], Ks))[m] for i in 1:Nint, m in 1:nwannier], "-k"
    )
    plot(
        xint .+ 2,
        [sort(λ([0, xint[i], xint[i]], Ks))[m] for i in 1:Nint, m in 1:nwannier],
        "-k",
    )

    #original
    plot(
        x,
        [sort(real(eigen(H[:, :, i, 1]).values))[m] for i in 1:N, m in 1:1],
        "o";
        label="Samples",
    )
    gca()[:set_color_cycle](nothing)
    plot(x, [sort(real(eigen(H[:, :, i, 1]).values))[m] for i in 1:N, m in 1:nwannier], "o")
    gca()[:set_color_cycle](nothing)
    plot(
        x .+ 1,
        [sort(real(eigen(H[:, :, 1, i]).values))[m] for i in 1:N, m in 1:nwannier],
        "o",
    )
    gca()[:set_color_cycle](nothing)
    plot(
        x .+ 2,
        [sort(real(eigen(H[:, :, i, i]).values))[m] for i in 1:N, m in 1:nwannier],
        "o",
    )
    gca()[:set_color_cycle](nothing)

    #interpolated
    plot(
        xint,
        [sort(real(eigen(Hint[:, :, i, 1]).values))[m] for i in 1:Nint, m in 1:1],
        "-";
        label="Interpolation",
    )
    gca()[:set_color_cycle](nothing)
    plot(
        xint,
        [sort(real(eigen(Hint[:, :, i, 1]).values))[m] for i in 1:Nint, m in 1:nwannier],
        "-",
    )
    gca()[:set_color_cycle](nothing)
    plot(
        xint .+ 1,
        [sort(real(eigen(Hint[:, :, 1, i]).values))[m] for i in 1:Nint, m in 1:nwannier],
        "-",
    )
    gca()[:set_color_cycle](nothing)
    plot(
        xint .+ 2,
        [sort(real(eigen(Hint[:, :, i, i]).values))[m] for i in 1:Nint, m in 1:nwannier],
        "-",
    )

    legend()

    xlabel("k")
    ylabel("ɛk")
    savefig("interp_2D.pdf")
end

if dim == 3
    # for iwann = 1:1
    #     figure()
    #     pl = pcolormesh(abs.(Wfour[iwann,div(end,4),:,:]),linewidth=0,rasterized=true) #force rasterization to avoid "kilt" effect when saving
    #     figure()
    #     pl = pcolormesh(abs.(Wfour[iwann,div(end,2),:,:]),linewidth=0,rasterized=true) #force rasterization to avoid "kilt" effect when saving
    #     figure()
    #     pl = pcolormesh(abs.(Wfour[iwann,div(3*end,4),:,:]),linewidth=0,rasterized=true) #force rasterization to avoid "kilt" effect when saving
    # end

    # xlim(-2,2)
    # ylim(-2,2)
    # xlabel("ξx")
    # ylabel("ξy")
    # savefig("2D_abs_$iwann.pdf",dpi=450)

    # data = abs((Wfour[1,:,:,:]))
    # using GLVisualize, GLWindow
    # window = glscreen()
    # timesignal = bounce(linspace(Float32(minimum(data)), Float32(maximum(data)),360))
    # vol = visualize(data, :iso, isovalue=timesignal)
    # _view(vol, window)
    # renderloop(window)

    iwann = 1
    arr = real(Wfour[iwann, :, :, :])
    grad = zeros(size(arr))
    for i in 1:(size(arr, 1) - 1)
        for j in 1:(size(arr, 2) - 1)
            for k in 1:(size(arr, 3) - 1)
                grad[i, j, k] = sqrt(
                    (arr[i + 1, j, k] - arr[i, j, k])^2 +
                    (arr[i, j + 1, k] - arr[i, j, k])^2 +
                    (arr[i, j, k + 1] - arr[i, j, k])^2,
                )
            end
        end
    end
    data = grad
    using GLVisualize, GLWindow
    window = glscreen()
    timesignal = bounce(linspace(Float32(minimum(data)), Float32(maximum(data)), 360))
    vol = visualize(data, :iso; isovalue=timesignal)
    _view(vol, window)
    renderloop(window)
end
