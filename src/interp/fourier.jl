using LinearAlgebra

## NEW
k_cryst(k) = k

"""
    fourier(f::Function, q_vectors, R_vectors)

Performs a fourier transform from the ab-initio kpoints to the wigner seitz unit cells.
The function will be executed inside the fourier transform loop, being called like
`f(iR, ik, phase)`
"""
function fourier(f::Function, q_vectors, R_vectors)
    for iR in 1:length(R_vectors)
        for ik in 1:length(q_vectors)
            phase = exp(-2im * π * (k_cryst(q_vectors[ik]) ⋅ R_vectors[iR]))
            f(iR, ik, phase)
        end
    end
end
    
#TODO Make a nice Fourier struct holding ik, iR, tbblock and factor
"Fourier transforms the tight binding hamiltonian and calls the R_function with the current index and the phase."
function invfourier(R_function::Function, tb_hami::TBHamiltonian{T}, kpoint::Vec3) where {T}
    
    for (iR, b) in enumerate(tb_hami)
        fac = ℯ^(2im * π * (b.R_cryst ⋅ kpoint))
        for i in eachindex(block(b))
            R_function(i, iR, b.R_cart, b, fac)
        end
    end
end

