using LinearAlgebra
using Dates


# Actualize the overlap matrix (if A is changed, the overlap should be changed accordingly)
function actualize_M!(p)
    M = p.M
    for K1 = 1:(p.N1*p.N2*p.N3)
        for nneighbor = 1:p.nntot
            i1, j1, k1 = p.K_to_ijk[K1, :]
            i2, j2, k2 = p.K_to_ijk[p.neighbors[K1, nneighbor], :]
            A1 = normalize_matrix(p.A[:, :, i1, j1, k1])
            A2 = normalize_matrix(p.A[:, :, i2, j2, k2])
            M[:, :, K1, nneighbor] = A1' * A2
        end
    end
    return M
end


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

function interpfft(A, newSize)
    Y = fft(A)
    (N1, N2, N3) = size(A)
    ##Assuming N1,N2,N3 are even
    Y_padded = zeros(ComplexF64, newSize, newSize, newSize)
    Y_padded[div(newSize - N1, 2)+1:div(newSize + N1, 2), div(newSize - N2, 2)+1:div(newSize + N2, 2), div(newSize - N3, 2)+1:div(newSize + N3, 2)] = fftshift(Y)
    Y_padded = ifftshift(Y_padded)
    Ainterp = (newSize^3 / (N1 * N2 * N3)) * real(ifft(Y_padded))
    return Ainterp
end

function interpFourier(A)
    Y = fftshift(fft(A))
    half_size = floor.(size(A) ./ 2)
    function interpA(x)
        poly = 0
        for k = 1:size(A)[3], j = 1:size(A)[2], i = 1:size(A)[1]
            poly += Y[i, j, k] * exp(im * 2 * pi * dot([(i - 1); (j - 1); (k - 1)] .- half_size, x))
        end
        return real(poly / length(A))
    end
    return interpA
end
