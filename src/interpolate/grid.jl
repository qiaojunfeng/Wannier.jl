using LinearAlgebra


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
