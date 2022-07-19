using LinearAlgebra

"""
Propagate A, defined at the first kpt, to the given list of kpts.
Those must be neighbors, and only the first kpoint is assumed to have been rotated.
"""
function propagate!(
    A::Array{T,3}, kpts::Vector{Int}, M::Array{T,4}, kpb_k::Matrix{Int}
) where {T<:Complex}
    N = length(kpts)
    n = size(A, 1)
    m = size(A, 2)
    @assert n == m "Non square matrix given as argument in propagate"

    for i in 2:N
        ik = kpts[i]
        ik0 = kpts[i - 1]
        Mᵏᵇ = overlap(M, kpb_k, ik, ik0)
        A[:, :, ik] = orthonorm_lowdin(Mᵏᵇ * A[:, :, ik0])
    end

    return nothing
end

"""
Choose the column and the target point to contract the
vector path of the columns of matrix_path
"""
function choose_pole(
    matrix_path::Array{T,3}, columns::Vector{Int}, prev_poles::Matrix{T}
) where {T<:Complex}
    # Number of iterations to find a pole far enough from the path
    n_iter = 100

    if isempty(columns)
        columns = [1]
    else
        columns = [columns; length(columns) + 1]
    end

    n_col = size(matrix_path, 2)
    n_k = size(matrix_path, 3)

    col = columns[end]
    @assert col <= n_col

    vec_path = matrix_path[:, col, :]

    # diameter
    diam = maximum([norm(vec_path[:, i] - vec_path[:, j]) for i in 1:n_k, j in 1:n_k])
    @debug "diamater = $diam"

    pole = zero(vec_path[:, 1])
    pole[col] = 1.0

    m = 1
    max_dist = 0.0
    max_pole = deepcopy(pole)

    if diam > 1.5
        # If the diameter of the path on the sphere is too big,
        # try a list of cardinal points as poles.
        while m < n_iter
            dist = minimum([norm(pole + vec_path[:, i]) for i in 1:n_k])

            if dist > max_dist
                max_pole = deepcopy(pole)
                max_dist = dist
            end

            if m <= n_col
                pole = zero(pole)
                pole[m] = 1
            elseif m > n_col && m <= 2 * n_col
                pole = zero(pole)
                pole[m - n_col] = -1
            else
                pole = randn(n_col) + randn(n_col) * im
            end

            P = Matrix{T}(I, n_col, n_col) - prev_poles * prev_poles'
            pole = P * pole
            pole = pole / norm(pole)

            m += 1
        end

        pole = max_pole

        # @info "Pole chosen" pole max_dist
        if m <= 2 * n_col
            @info "Pole chosen by cardinal point" max_dist
        else
            @info "Pole chosen randomly" max_dist
        end
    else
        # If the diameter of the path is small, the Barycenter of the path is a good pole.
        @info "Pole chosen by Barycenter of path"
        pole = dropdims(sum(vec_path; dims=2); dims=2)

        @assert norm(pole) > 1e-1
        pole = pole / norm(pole)
    end

    return columns, pole
end

"""
Compute the parallel transport of the matrix path (matrix_path)
along the homotopy of the frame (frame_path) that starts with
the columns (columns[]) of matrix_path.

size(frame_path) = (n_col, n_col, n_k, n_t)
size(matrix_path) = (n_col, n_col, n_k)
"""
function matrix_parallel_transport(
    frame_path::Array{T,4},
    matrix_path::Array{T,3},
    columns::Vector{Int},
    backwards::Bool=false,
) where {T<:Complex}
    n_k = size(frame_path, 3)
    n_t = size(frame_path, 4)
    n_col = size(matrix_path, 2)

    @assert size(frame_path, 1) == n_col
    @assert size(matrix_path) == (n_col, n_col, n_k)

    P = zeros(T, n_col, n_col, n_k, n_t)
    U = zeros(T, n_col, n_col, n_k, n_t)

    # for ik = 1:n_k, ic = 1:n_col, it = 1:n_t
    # The final spread does not change much w.r.t. the order,
    # however the numbers in AMN file will be a little bit different.
    # I choose ic increases the fastest because it is more friendly for cache.
    for it in 1:n_t, ik in 1:n_k, ic in 1:n_col
        if !backwards
            # Compute parallel transport starting from it = 1 to it = n_t
            P[:, :, ik, it] = I - frame_path[:, :, ik, it] * frame_path[:, :, ik, it]'

            # Propagate column ic of matrix_path along frame_path
            if ic ∈ columns
                idx = findfirst(column -> column == ic, columns)

                if abs(norm(frame_path[:, idx, ik, it]) - 1) < 1e-2
                    U[:, ic, ik, it] = frame_path[:, idx, ik, it]
                else
                    if it > 1
                        U[:, ic, ik, it] = P[:, :, ik, it] * U[:, ic, ik, it - 1]
                    else
                        U[:, ic, ik, it] = P[:, :, ik, it] * matrix_path[:, ic, ik]
                    end
                end
            else
                if it > 1
                    U[:, ic, ik, it] = P[:, :, ik, it] * U[:, ic, ik, it - 1]
                else
                    U[:, ic, ik, it] = P[:, :, ik, it] * matrix_path[:, ic, ik]
                end
            end

            U[:, ic, ik, it] /= norm(U[:, ic, ik, it])
        else
            # Compute parallel transport from it = n_t to it = 1
            iit = n_t - it + 1  # inverse of it
            P[:, :, ik, iit] = I - frame_path[:, :, ik, iit] * frame_path[:, :, ik, iit]'

            # Propagate column ic of matrix_path along frame_path
            if ic ∈ columns
                idx = findfirst(column -> column == ic, columns)

                if abs(norm(frame_path[:, idx, ik, iit]) - 1) < 1e-2
                    U[:, ic, ik, iit] = frame_path[:, idx, ik, iit]
                else
                    if it > 1
                        U[:, ic, ik, iit] = P[:, :, ik, iit] * U[:, ic, ik, iit + 1]
                    else
                        U[:, ic, ik, iit] = P[:, :, ik, iit] * matrix_path[:, ic, ik]
                    end
                end
            else
                if it > 1
                    U[:, ic, ik, iit] = P[:, :, ik, iit] * U[:, ic, ik, iit + 1]
                else
                    U[:, ic, ik, iit] = P[:, :, ik, iit] * matrix_path[:, ic, ik]
                end
            end

            U[:, ic, ik, iit] /= norm(U[:, ic, ik, iit])
        end

        not_columns = [j for j = 1:n_col if j ∉ columns]

        # Normalize the remaining columns
        if length(not_columns) > 0
            U[:, not_columns, ik, it] = orthonorm_lowdin(U[:, not_columns, ik, it])
        end
    end

    return U
end

"""
Create an interpolation path between x and y, unit vectors
by normalizing the linear interpolation of parameter t[i].
"""
function interpolate_vec(
    x::AbstractVector{T}, y::AbstractVector{T}, t::AbstractVector{FT}
) where {T<:Union{Complex,Real},FT<:Real}
    @assert length(t) > 2
    @assert length(x) == length(y)

    n_x = length(x)
    n_t = length(t)

    v = zeros(T, n_x, n_t)

    for i in 1:n_t
        v[:, i] = (1 - t[i]) * x + t[i] * y

        n = norm(v[:, i])
        @assert n > 1e-2

        v[:, i] /= n
    end

    return v
end

"""
Contract the matrix_path to a single matrix point
in the space of unitaries using parallel transport.

i.e., contract Obs(k) matrices to constant vectors.

size(matrix_path) = n_wann x n_wann x n_k
t: vector of kpoint indexes along a different k direction
"""
function matrix_transport(matrix_path::Array{Complex{T},3}, t::Vector{T}) where {T<:Real}
    n_row, n_col, n_k = size(matrix_path)
    @assert n_row == n_col

    # interpolate matrix along t direction
    n_t = length(t)

    if n_t == 1
        # the interpolated matrices are stored along the 4th dimension
        return reshape(matrix_path, (n_row, n_col, n_k, 1))
    end

    @assert begin
        for i in 1:n_k
            P = matrix_path[:, :, i] * matrix_path[:, :, i]'
            norm(P - I) >= 1e-5 && return false
        end
        true
    end

    CT = Complex{T}
    columns = zeros(Int, 0)
    prev_poles = zeros(CT, n_col, n_col)
    # frame along the path
    frame_path = zeros(CT, n_col, n_col, n_k, n_t)
    U = zeros(CT, n_col, n_col, n_k, n_t)

    for col in 1:n_col
        # TODO ensure pole' * U[:, col+1, i, end] != -1
        columns, pole = choose_pole(matrix_path, columns, prev_poles)
        prev_poles[:, col] = pole
        @debug "matrix_transport" columns pole

        U = matrix_parallel_transport(frame_path, matrix_path, columns, true)

        # Compute the coefficients of the transported path on the basis of the poles,
        # and contract it to a point.
        if col < n_col
            n_c = n_col - col + 1
            c = zeros(CT, n_c, n_k, n_t)

            e1 = zeros(CT, n_c)
            e1[1] = 1

            for ik in 1:n_k
                Us = zeros(CT, n_col, n_c, n_k)

                # compute GLS2019 Eq. (7): c(k,t=1)
                for j in 1:n_c
                    Us[:, j, ik] = U[:, col + j - 1, ik, end]
                    c[j, ik, 1] = Us[:, j, ik]' * pole
                end

                # from c(k,t=1) to c(k,t)
                for j in 1:n_t
                    c[:, ik, j] = (1 - t[j]) * c[:, ik, 1] + t[j] * e1
                    c[:, ik, j] /= norm(c[:, ik, j])
                end

                # new frame_path
                for j in 1:n_t
                    frame_path[:, col, ik, j] = sum(
                        U[:, col + k - 1, ik, j] * c[k, ik, j] for k in 1:n_c
                    )
                    frame_path[:, col, ik, j] /= norm(frame_path[:, col, ik, j])
                end
            end
        else
            # Contract the last column by finding a continuous phase.
            ϕ = zeros(T, n_k)

            for ik in 1:n_k
                ϕ[ik] = imag(log(pole' * U[:, col, ik, 1]))

                if ik > 1
                    kmin = argmin([abs(ϕ[ik] + 2π * k - ϕ[ik - 1]) for k in -1:1])
                    ϕ[ik] += (kmin - 2) * 2π
                end
            end

            m = (ϕ[end] - ϕ[1]) / 2π
            # TODO
            # @assert m % 1 ≈ 0 "m = $m"
            @info "Chern number" m

            for i in 1:n_k, j in 1:n_t
                U[:, col, i, j] *= exp(-im * (1 - t[j]) * ϕ[i])
            end
        end
    end

    # Bring the contraction point from Obs to I
    O = orthonorm_lowdin(U[:, :, 1, 1])
    @debug "obstruction matrix" O

    for i in 1:n_k
        for j in 1:n_t
            U[:, :, i, j] = powm(O', 1 - t[j]) * U[:, :, i, j]
            # U[:, :, i, j] = orthonorm_lowdin(U[:, :, i, j])
        end
    end

    return U
end

struct Obstruction{T<:Complex}
    # obstruction matrix in x-y plane, at ky = 1 along kx = 0 -> 1
    Oxy::Array{T,3}

    # obstruction matrix in x-z plane, at kz = 1 along kx = 0 -> 1
    Oxz::Array{T,3}

    # obstruction matrix in y-z plane, at kz = 1 along ky = 0 -> 1
    Oyz::Array{T,3}

    # frame in x-y plane, at ky = 1 along kx = 0 -> 1
    Uxy::Array{T,4}

    # frame in x-z plane, at kz = 1 along kx = 0 -> 1
    Uxz::Array{T,4}

    # frame in y-z plane, at kz = 1 along ky = 0 -> 1
    Uyz::Array{T,4}
end
