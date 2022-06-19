using LinearAlgebra

#Create an interpolation path between x1 and x2, unit vectors
#by normalizing the linear interpolation of parameter t[i].
function interpolate_vec(x1, x2, t)
    @assert(size(t)[1] > 2)
    v_interp = zeros(ComplexF64, size(x1)[1], size(t)[1])
    for i in 1:size(t)[1]
        v_interp[:, i] = (1 - t[i]) * x1 + t[i] * x2
        n = norm(v_interp[:, i])
        @assert(n > 1e-2)
        v_interp[:, i] = v_interp[:, i] / n
    end
    return v_interp
end


#Compute the parallel transport of the matrix path matrixPath
#along the homotopy of the frame framePath that starts with
#the columns columns[] of matrixPath.
function parallelTransport(framePath, columns, matrixPath, backwards=false)
    @assert(size(matrixPath)[1] == size(framePath)[1])
    nc = size(columns)[1]
    ntot = size(matrixPath)[1]
    n = size(matrixPath)[end]
    P = zeros(ComplexF64, ntot, ntot, n, n)
    U = zeros(ComplexF64, ntot, ntot, n, n)
    @assert(size(framePath[1, 1, :, :] * framePath[1, 1, :, :]') == (n, n))

    for i1 = 1:ntot, i2 = 1:ntot, j = 1:n
        # Compute parallel transport starting from i2=1 to i2=ntot
        if !backwards
            P[i1, i2, :, :] = I - framePath[i1, i2, :, :] * framePath[i1, i2, :, :]'
            # Propagate column j of matrixPath along framePath
            if j in columns
                ind = find(column -> column == j, columns)
                if (abs(norm(framePath[i1, i2, :, ind]) - 1) < 1e-2)
                    U[i1, i2, :, j] = framePath[i1, i2, :, ind]
                else
                    if i2 > 1
                        U[i1, i2, :, j] = P[i1, i2, :, :] * U[i1, i2-1, :, j]
                    else
                        U[i1, i2, :, j] = P[i1, i2, :, :] * matrixPath[i1, :, j]
                    end
                end
            else
                if i2 > 1
                    U[i1, i2, :, j] = P[i1, i2, :, :] * U[i1, i2-1, :, j]
                else
                    U[i1, i2, :, j] = P[i1, i2, :, :] * matrixPath[i1, :, j]
                end
            end
            U[i1, i2, :, j] /= norm(U[i1, i2, :, j])
            # Compute parallel transport from i2=ntot to i2=1
        else
            P[i1, ntot-i2+1, :, :] = I - framePath[i1, ntot-i2+1, :, :] * framePath[i1, ntot-i2+1, :, :]'
            # Propagate column j of matrixPath along framePath
            if j in columns
                ind = (findall(column -> column == j, columns))[1]
                if (abs(norm(framePath[i1, ntot-i2+1, :, ind]) - 1) < 1e-2)
                    U[i1, ntot-i2+1, :, j] = framePath[i1, ntot-i2+1, :, ind]
                else
                    if i2 > 1
                        U[i1, ntot-i2+1, :, j] = P[i1, ntot-i2+1, :, :] * U[i1, ntot-i2+2, :, j]
                    else
                        U[i1, ntot-i2+1, :, j] = P[i1, ntot-i2+1, :, :] * matrixPath[i1, :, j]
                    end
                end
            else
                if i2 > 1
                    U[i1, ntot-i2+1, :, j] = P[i1, ntot-i2+1, :, :] * U[i1, ntot-i2+2, :, j]
                else
                    U[i1, ntot-i2+1, :, j] = P[i1, ntot-i2+1, :, :] * matrixPath[i1, :, j]
                end
            end
            U[i1, ntot-i2+1, :, j] /= norm(U[i1, ntot-i2+1, :, j])

        end
        notColumns = [!(j in columns) for j = 1:n]
        # Normalize the remaining columns
        if (!all(.!notColumns))
            U[i1, i2, :, notColumns] = normalize_matrix(U[i1, i2, :, notColumns])
        end

    end
    return U
end


#Choose the column and the target point to contract the
#vector path of the columns of matrixPath
function choosePole(columns, matrixPath, t, previousPoles)
    #Number of iterations to find a pole far enough from the path
    nbIt = 100
    if isempty(columns)
        columns = [1]
    else
        columns = [columns; (length(columns) + 1)]
    end
    col = columns[end]
    @assert(col <= size(matrixPath)[2])
    n = size(matrixPath)[2]
    vecPath = matrixPath[:, col, :]
    diam = maximum([norm(vecPath[i, :] - vecPath[j, :]) for i = 1:size(t)[1], j = 1:size(t)[1]])
    pole = deepcopy(vecPath[1, :] * 0.0)
    pole[col] = 1.0
    m = 1
    maxDist = 0.0
    maxPole = pole

    # If the diameter of the path on the sphere is too big,
    # try a list of cardinal points as poles.
    if (diam > 1.5)
        while (m < nbIt)
            dist = minimum([norm(pole + vecPath[i, :]) for i = 1:size(t)[1]])

            if dist > maxDist
                maxPole = pole
                maxDist = dist
            end

            if (m <= n)
                pole = Matrix((1.0 + 0.0im)I, n, n)[m, :]
            elseif (m > n && m <= 2n)
                pole = -Matrix((1.0 + 0.0im)I, n, n)[m-n, :]
            else
                pole = [randn() + randn() * im for i = 1:size(pole)[1]]
            end
            Proj = Matrix((1.0 + 0.0im)I, length(pole), length(pole)) - previousPoles * previousPoles'
            pole = Proj * pole
            pole = pole / norm(pole)

            m += 1
        end
        pole = maxPole
        # If the diameter of the path is small, the barycentre of the path is a good pole.
    else
        println("Pole chosen by barycentre of path")
        pole = sum(vecPath[i, :] for i = 1:size(t)[1])
        @assert(norm(pole) > 1e-1)
        pole = pole / norm(pole)
    end

    #println("Pole chosen = $pole\nDistance to nearest point = $maxDist")
    return columns, pole

end


#Contract the matrix path matrixPath to a single matrix point
#in the space of unitaries using parallel transport
function matrixTransport(matrixPath, t)
    if (size(t)[1] == 1)
        return reshape(matrixPath, (size(matrixPath)[1], 1, size(matrixPath)[2], size(matrixPath)[3]))
    end

    @assert(maximum([norm(matrixPath[i, :, :] * matrixPath[i, :, :]' - I) for i = 1:size(t)[1]]) < 1e-5)
    n = size(matrixPath)[end]
    columns = zeros(Int64, 0)
    previousPoles = zeros(ComplexF64, n)
    vecInterp = zeros(ComplexF64, size(t)[1], size(t)[1], n, n)
    U = zeros(ComplexF64, size(t)[1], size(t)[1], n, n)
    for col = 1:n
        ### ensure pole'*U[i,end,:,col+1] != -1
        columns, pole = choosePole(columns, matrixPath, t, previousPoles)
        println("columns, pole = $columns,  $pole")
        #vecPath = matrixPath[:,:,col]

        U = parallelTransport(vecInterp, columns, matrixPath, true)
        # Compute the coefficients of the transported path on the basis of the poles,
        # and contract it to a point.
        if col < n
            c = zeros(ComplexF64, size(t)[1], size(t)[1], n - col + 1)
            notColumns = [!(j in columns[1:end-1]) for j = 1:n]
            for i = 1:size(t)[1]
                Us = zeros(ComplexF64, size(t)[1], n - col + 1, n)
                for k = 1:size(c)[end]
                    Us[i, k, :] = U[i, end, :, col+k-1]
                    c[i, 1, k] = Us[i, k, :]' * pole
                end
                for j = 1:size(t)[1]
                    e1 = Matrix((1.0 + 0.0im)I, n - col + 1, n - col + 1)[1, :]
                    c[i, j, :] = (1 - t[j]) * c[i, 1, :] + t[j] * e1
                    c[i, j, :] = c[i, j, :] / norm(c[i, j, :])
                end
                for j = 1:size(t)[1]
                    vecInterp[i, j, :, col] = sum(U[i, j, :, col+k-1] * c[i, j, k] for k = 1:size(c)[end])
                    vecInterp[i, j, :, col] /= norm(vecInterp[i, j, :, col])
                end
            end
            # Contract the last column by finding a continuous phase.
        else
            logPhase = zeros(Float64, size(t)[1])
            for i = 1:size(t)[1]
                #println("U[$i,1,:,$(col-1)] = $(U[i,1,:,col-1])")
                #println("U[$i,1,:,$col] = $(U[i,1,:,col])")
                logPhase[i] = imag(log(pole' * U[i, 1, :, col]))
                if i > 1
                    kmin = argmin([abs(logPhase[i] + 2 * pi * k - logPhase[i-1]) for k in -1:1])
                    logPhase[i] += (kmin - 2) * 2 * pi
                end
            end
            for i = 1:size(t)[1], j = 1:size(t)[1]
                U[i, j, :, col] = U[i, j, :, col] * exp(-im * (1 - t[j]) * logPhase[i])
            end
        end
    end
    # Bring the contraction point from Obs to I
    Obs = normalize_matrix(U[1, 1, :, :])
    println("Obs = $Obs")
    for i = 1:size(t)[1]
        for j = 1:size(t)[1]
            U[i, j, :, :] = powm(Obs', 1 - t[j]) * U[i, j, :, :]
            #U[i,j,:,:] = normalize_matrix(U[i,j,:,:])
        end
    end
    return U
end
