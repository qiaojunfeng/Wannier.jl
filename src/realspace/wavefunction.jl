
using Printf: @sprintf
using LazyGrids: ndgrid

export read_realspace_wf, write_realspace_wf

"""
    _get_unk_ext(unkdir::AbstractString)

Get the extension name of `UNK` files.

For example:
- `1` of `UNK00001.1`, for no-spin calcluation
- `NC` of `UNK00001.NC`, for non-collinear calculation
"""
function _get_unk_ext(unkdir::AbstractString)
    files = readdir(unkdir)
    i = findfirst(x -> startswith(x, "UNK"), files)
    isnothing(i) && error("No UNK files found in $unkdir")
    suffix = split(files[i], ".")[2]  # "1" or "NC"
    return suffix
end

"""
    read_realspace_wf(U, kpoints, n_supercells, unkdir="."; R=[0, 0, 0])

Read `UNK` files, rotate gauge, and generate real space WFs.

# Arguments
- `U`: `n_bands * n_wann * n_kpts`, gauge rotation matrix
- `kpoints`: `3 * n_kpts`, each column is a vector for fractional coordinates
- `n_supercells`: an integer or a vector of 3 integers for 3 directions along lattice vectors,
    defines the number of super cells where the real space WF lives
- `unkdir`: folder containing `UNK` files

# Keyword arguments
- `R`: fractional coordinates w.r.t lattice (actually integers), the cell for WF
    ``|n \\bm{R} \\rangle``, default is generating WF at home unit cell ``|n \\bm{0} \\rangle``

!!! tip

    See also the section [Normalization convention of WFs](@ref) for further explanation.
"""
function read_realspace_wf(
    U::Vector{Matrix{Complex{T}}},
    kpoints::Vector{Vec3{T}},
    lattice,
    n_supercells::AbstractVector{Int},
    unkdir::AbstractString;
    R::AbstractVector{Int}=[0, 0, 0],
    wan_plot_list = 1:size(U[1], 2)
) where {T<:Real}
    n_bands, n_wann = size(U[1])
    
    nwfun = length(wan_plot_list)
    
    n_kpts = length(U)
    length(n_supercells) == 3 || error("n_supercells must be 3-vector")

    supercells = Vector{UnitRange}()
    for i in 1:3
        n_s = n_supercells[i]
        d, r = divrem(n_s, 2)
        if r == 0
            push!(supercells, (-d):(d - 1))
        else
            push!(supercells, (-d):d)
        end
    end
    n_sx, n_sy, n_sz = n_supercells

    ext = _get_unk_ext(unkdir)
    unkname(ik) = joinpath(unkdir, @sprintf("UNK%05d.%s", ik, ext))

    # read 1st kpoint to allocate matrices
    ik = 1
    _, Ψₖ = read_unk(unkname(ik))

    n_gx, n_gy, n_gz, _, ns = size(Ψₖ)
    size(Ψₖ, 4) == n_bands || error("incompatible n_bands")

    # WF in realspace
    W = zeros(eltype(U[1]), n_gx * n_sx, n_gy * n_sy, n_gz * n_sz, ns, nwfun)

    # generate X, Y, Z fractional coordinates relative to lattice (here unknown)
    # actually X./n_gx, Y./n_gy, Z./n_gz are the fractional coordinates w.r.t lattice
    X = (supercells[1][1] * n_gx):((supercells[1][end] + 1) * n_gx - 1)
    Y = (supercells[2][1] * n_gy):((supercells[2][end] + 1) * n_gy - 1)
    Z = (supercells[3][1] * n_gz):((supercells[3][end] + 1) * n_gz - 1)
    # we start from origin, however W90 starts from the left of origin
    # uncomment these lines to exactly reproduce W90 output
    # seems subtracting 1 leads to much worse performance
    #  11.015472 seconds (136.09 M allocations: 5.570 GiB, 10.77% gc time)
    #   3.513259 seconds (25.60 M allocations: 2.728 GiB, 24.68% gc time)
    # X = X .- 1
    # Y = Y .- 1
    # Z = Z .- 1

    """Modify W"""
    @inbounds function add_k!(ik, Ψₖ)
        k = kpoints[ik]
        # rotate Ψ
        # * does not support high dimensional matrix multiplication,
        # I need to reshape it to 2D matrix
        for is = 1:ns
            ΨUₖ = reshape(reshape(view(Ψₖ, :, :, :, :, is), :, n_bands) * view(U[ik], :, wan_plot_list), n_gx, n_gy, n_gz, nwfun)
            for iw in 1:nwfun
                # make sure Ψ is normalized to 1 in unit cell
                # the QE output UNK needs this factor to be normalized
                # Note in QE, the r->G FFT has a factor 1/N,
                # while the G->r IFFT has factor 1,
                # but FFT/IFFT is unitary only if they have factor sqrt(1/N)
                Threads.@threads for iz in 1:length(Z)
                    z = Z[iz]
                    # ΨUₖ is only defined in the home unit cell, find corresponding indexes
                    # e.g. x=0 -> W[1,1,1], x=1 -> W[2,1,1], x=n_gx -> W[1,1,1], x=-1 -> W[end,1,1]
                    jz = mod(z, n_gz) + 1
                    for (iy, y) in enumerate(Y)
                        jy = mod(y, n_gy) + 1
                        for (ix, x) in enumerate(X)
                            jx = mod(x, n_gx) + 1
                            # r grid in fractional coordinates
                            r = Vec3(x / n_gx, y / n_gy, z / n_gz)
                            # UNK file (and ΨUₖ) is the periodic part of Bloch wavefunction,
                            # need a factor exp(ikr)
                            f = exp(2π * im * (k ⋅ (r - R)))
                            W[ix, iy, iz, is, iw] += f * ΨUₖ[jx, jy, jz, iw]
                        end
                    end
                end
            end
        end
        return nothing
    end

    # 1st kpoint
    add_k!(ik, Ψₖ)

    for ik in 2:n_kpts
        unk = unkname(ik)
        _, Ψₖ = read_unk(unk)
        add_k!(ik, Ψₖ)
    end

    W ./= (n_kpts * sqrt(n_gx * n_gy * n_gz))
    
    if ns == 1
        for i in 1:nwfun
            f = fix_global_phase(W[:, :, :, 1, i])
            @views W[:, :, :, 1, i] .*= f
        end
    end
    
    for i in 1:nwfun
        for is in 1:ns
            r = compute_imre_ratio(W[:, :, :, is, i])
            @printf("Max Im/Re ratio of WF %4d; spin %4d = %10.4f\n", i, is, r)
        end
    end

    points = [
        lattice * Vec3((x - 1) / n_sx, (y - 1) / n_sy, (z - 1) / n_sz) for
        x in X, y in Y, z in Z
    ]
   
    if ns == 1
        wfuncs_out = Vector{WannierFunction{1,eltype(W).parameters[1]}}(
            undef, nwfun
        )
        Threads.@threads for i in 1:nwfun
            wfuncs_out[i] = WannierFunction{1,eltype(W).parameters[1]}(
                points, map(x -> SVector(x), view(W, :, :, :, 1, i))
            )
        end
    else
        wfuncs_out = Vector{WannierFunction{2,eltype(W).parameters[1]}}(
            undef, nwfun
        )
        Threads.@threads for i in 1:size(wfuncs_all, 1)
            wfuncs_out[i] = WannierFunction{2,eltype(W).parameters[1]}(
                points,
                map(
                    x -> SVector(x),
                    zip(view(W, :, :, :, 1, i), view(W, :, :, :, 2, i)),
                ),
            )
        end
    end
    # to fractional coordinates
    Xg, Yg, Zg = ndgrid(X ./ n_gx, Y ./ n_gy, Z ./ n_gz)
    return Xg, Yg, Zg, wfuncs_out
end

function read_realspace_wf(
    U::Vector{Matrix{Complex{T}}},
    kpoints::Vector{Vec3{T}},
    lattice,
    n_supercells::Int,
    unkdir::AbstractString;
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    return read_realspace_wf(
        U, kpoints, lattice, [n_supercells, n_supercells, n_supercells], unkdir; R=R
    )
end

"""
    read_realspace_wf(lattice, U, kpoints, n_supercells=2, unkdir="."; R=[0, 0, 0])

Read `UNK` files, rotate gauge, and generate real space WFs.

This is a more user-friendly version, which returns a tuple of `(RGrid, W)`,
where `RGrid` is the grid on which `W` is defined, and `W` is volumetric data for WFs.

# Arguments
- `lattice`: each column is a lattice vector
- `U`: `n_bands * n_wann * n_kpts`, gauge rotation matrix
"""
function read_realspace_wf(
    lattice::AbstractMatrix{T},
    U::Vector{Matrix{Complex{T}}},
    kpoints::Vector{Vec3{T}},
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    X, Y, Z, W = read_realspace_wf(U, kpoints, lattice, n_supercells, unkdir; R=R)
    rgrid = RGrid(lattice, X, Y, Z)
    return rgrid, W
end

"""
    read_realspace_wf(model, U, n_supercells=2, unkdir="."; R=[0, 0, 0])

Read `UNK` files, rotate gauge, and generate real space WFs.

This is a most user-friendly version, use `lattice`, `U` and `kpoints` from `model`
and returns a tuple of `(RGrid, W)`,
where `RGrid` is the grid on which `W` is defined, and `W` is volumetric data for WFs.

# Arguments
- `model`: a `Model`
- `U`: `n_bands * n_wann * n_kpts`, gauge rotation matrix
"""
function read_realspace_wf(
    model::Model{T},
    U::Vector{Matrix{Complex{T}}},
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    return read_realspace_wf(model.lattice, U, model.kpoints, n_supercells, unkdir; R=R)
end

"""
    read_realspace_wf(model, n_supercells=2, unkdir="."; R=[0, 0, 0])

Read `UNK` files, rotate gauge, and generate real space WFs.

This is a most user-friendly version, use `lattice`, `U` and `kpoints` from `model`
and returns a tuple of `(RGrid, W)`,
where `RGrid` is the grid on which `W` is defined, and `W` is volumetric data for WFs.

# Arguments
- `model`: a `Model`
"""
function read_realspace_wf(
    model::Model{T},
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".";
    R::AbstractVector{Int}=[0, 0, 0],
) where {T<:Real}
    return read_realspace_wf(model, model.U, n_supercells, unkdir; R=R)
end

"""
    write_realspace_wf(seedname, U, kpoints, lattice, atom_positions, atom_labels;
        n_supercells=2, unkdir=".", part=real, format=:xsf, wf_center=nothing)

Write real space WFs to `xsf` or `cube` files.

# Arguments
- `seedname`: the name prefix for `cube` files, e.g., `seedname_00001.cube`
- `U`: gauge rotation matrix
- `kpoints`: each column is a kpoint, fractional coordinates
- `lattice`: each column is a lattice vector
- `atom_positions`: each column is an atom position, fractional coordinates w.r.t. lattice
- `atom_labels`: each element is an atom label

# Keyword arguments
- `n_supercells`: number of supercells in each direction,
    equivalent to `wannier_plot_supercell` of `Wannier90`
- `unkdir`: directory of `UNK` files
- `part`: which part to plot? pass a `Function`, e.g. `real`, `imag`, `abs2`
- `format`: `:xsf` or `:cube`
- `wf_center`: WF centers in fractional coordinates w.r.t. lattice.
    Only used for `cube` format, add additional atoms around WF centers.

!!! note

    `Wannier90` only plot the real part of WFs, so `part=real` is the default.

!!! tip

    See also the section [Normalization convention of WFs](@ref) for further explanation.
"""
function write_realspace_wf(
    seedname::AbstractString,
    U::AbstractArray,
    kpoints::Vector,
    lattice::AbstractMatrix,
    atom_positions::AbstractVector,
    atom_labels::AbstractVector{String};
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".",
    part::Function=real,
    format::Symbol=:xsf,
    wf_center::Union{Nothing,AbstractMatrix}=nothing,
)
    format ∈ [:xsf, :cube] || error("format must be :xsf or :cube")

    rgrid, W = read_realspace_wf(lattice, U, kpoints, n_supercells, unkdir)
    n_wann = length(W)
    ns = nspin(W[1])
    if format == :cube && wf_center === nothing
        # if not given, compute in realspace, lower accuracy. in fractional coordinates
        wf_center = map(c -> inv(lattice) * c, center(rgrid, W))
    end

    # In principle, since we normalize Bloch wavefunction inside one unit cell,
    # thus the WF is normalized inside the n_kpts times unit cell.
    # However, W90 does not check the normalization of UNK files, it just plainly
    # mulitply the unitary matrices. But in read_readspace_wf, we normalized
    # the unk files, so here we need to multiply a factor to reproduce the
    # W90 output XSF or cube files.
    if isa(n_supercells, Integer)
        n_supcells = n_supercells^3
    else
        n_supcells = prod(n_supercells)
    end
    W .*= sqrt(length(W[1].values) / n_supcells)


    # seems W90 always write the real part, so I use real as default
    atom_numbers = get_atom_number(atom_labels)

    if format == :xsf
        for i in 1:n_wann
            for is = 1:ns
                filename = ns == 1 ? @sprintf("%s_%05d.%s", seedname, i, "xsf") : @sprintf("%s_%05d_%05d.%s", seedname, i,is, "xsf")
                write_xsf(
                    filename, lattice, atom_positions, atom_numbers, rgrid, map(x -> part(x[1]), W[i].values)
                )
            end
        end
    elseif format == :cube
        for i in 1:n_wann
            for is = 1:ns
                filename = ns == 1 ? @sprintf("%s_%05d.%s", seedname, i, "cube") : @sprintf("%s_%05d_%05d.%s", seedname, i,is, "cube")
                write_cube(
                    filename,
                    lattice,
                    atom_positions,
                    atom_numbers,
                    wf_center[i],
                    rgrid,
                    map(x -> part(x[is]), W[i].values),
                )
            end
        end
    end

    return nothing
end

"""
    write_realspace_wf(seedname, model; n_supercells=2, unkdir=".", part=real, format=:xsf)

Write real space WFs to `xsf` or `cube` files.

This is a user-friendly version that use `model` to fill the arguments of `write_realspace_wf`.
"""
function write_realspace_wf(
    seedname::String,
    model::Model;
    n_supercells::Union{AbstractArray{Int},Int}=2,
    unkdir::AbstractString=".",
    part::Function=real,
    format::Symbol=:xsf,
)
    wf_center = nothing
    if format == :cube
        # compute in recip space, more accurate than realspace
        wf_center = map(c -> inv(model.lattice) * c, center(model))
    end
    return write_realspace_wf(
        seedname,
        model.U,
        model.kpoints,
        model.lattice,
        model.atom_positions,
        model.atom_labels;
        n_supercells,
        unkdir,
        part,
        format,
        wf_center,
    )
end

"""
Real space WannierFunction defined on a uniformly (in crystal coordinates) spaced r-grid. 
"""
struct WannierFunction{N,T<:AbstractFloat} <: AbstractArray{SVector{N,Complex{T}},3}
    points::Array{Vec3{T},3}
    values::Array{SVector{N,Complex{T}},3}
end

function WannierFunction(point_func::Function, points::Array)
    return WannierFunction(points, point_func.(points))
end

for f in (:size, :getindex, :setindex!)
    @eval @inline Base.@propagate_inbounds function Base.$f(x::WannierFunction, i...)
        return Base.$f(x.values, i...)
    end
end

for f in (:length, :stride, :ndims, :axes, :strides)
    @eval @inline Base.$f(w::WannierFunction) = Base.$f(w.values)
end

function Base.similar(x::WannierFunction, ::Type{S}) where {S}
    return WannierFunction(x.points, similar(x.values, S))
end

Base.unsafe_convert(T::Type{<:Ptr}, x::WannierFunction) = Base.unsafe_convert(T, x.values)

Base.Broadcast.broadcastable(w::WannierFunction) = w.values

nspin(w::WannierFunction{N}) where {N} = N

#### LinearAlgebra overloads
function LinearAlgebra.adjoint(w::WannierFunction)
    out = WannierFunction(w.points, similar(w.values))
    return adjoint!(out, w)
end

LinearAlgebra.adjoint!(w1::WannierFunction, w2::WannierFunction) = w1 .= adjoint.(w2)

function LinearAlgebra.dot(w1::WannierFunction{T}, w2::WannierFunction{T}) where {T}
    s = zero(T)
    for (v1, v2) in zip(w1.values, w2.values)
        s += v1' * v2
    end
    return real(s)
end

function LinearAlgebra.dot(v::Vector, wfs::Vector{<:WannierFunction})
    res = WannierFunction(wfs[1].points, zeros(eltype(wfs[1].values), size(wfs[1].values)))
    for ic in 1:length(v)
        res .+= v[ic] .* wfs[ic]
    end
    return res
end
LinearAlgebra.dot(wfs::Vector{<:WannierFunction}, v::Vector) = dot(v, wfs)

LinearAlgebra.norm(wfc::WannierFunction) = dot(wfc, wfc)

LinearAlgebra.normalize!(wfc::WannierFunction) = wfc ./= sqrt(norm(wfc))

same_grid(w1::WannierFunction, w2::WannierFunction) = w1.points === w2.points

function wan_op(op::Function, w1::W, w2::W) where {W<:WannierFunction}
    @assert same_grid(w1, w2) "Wannier functions are not defined on the same grid"
    return WannierFunction(w1.points, op(w1.values, w2.values))
end

Base.:(+)(w1::WannierFunction, w2::WannierFunction) = wan_op(+, w1, w2)
Base.:(*)(w1::WannierFunction, w2::WannierFunction) = wan_op(*, w1, w2)
Base.:(-)(w1::WannierFunction, w2::WannierFunction) = wan_op(-, w1, w2)
Base.:(*)(w1::WannierFunction, n::Number) = WannierFunction(w1.points, w1.values .* n)
Base.:(*)(n::Number, w1::WannierFunction) = WannierFunction(w1.points, n .* w1.values)
Base.:(/)(w1::WannierFunction, n::Number) = WannierFunction(w1.points, n ./ w1.values)
Base.:(/)(n::Number, w1::WannierFunction) = WannierFunction(w1.points, w1.values ./ n)

LinearAlgebra.dot(w1::WannierFunction, n::Number) = w1 * n
LinearAlgebra.dot(n::Number, w1::WannierFunction) = n * w1

function generate_wannierfunctions(
    model,
    unkdir::AbstractString,
    wannier_plot_supercell::NTuple{3,Int}=(3,3,3),
    wan_plot_list=1:size(U[1],2);
    R::Vec3{Int}=Vec3(0, 0, 0),
) where {T}
    return read_realspace_wf(model.U, model.kpoints, model.lattice, wannier_plot_supercell, unkdir; R, wan_plot_list)[2]
end

function bloch_sum(wfunc, kpoint; i_pos_offset=(0, 0, 0), i_neg_offset=(0, 0, 0))
    cell_boundaries = div.(size(wfunc.points), 3) .+ 1
    x = wfunc.points[cell_boundaries[1] + 1, 1, 1] .- wfunc.points[1]
    y = wfunc.points[1, cell_boundaries[2] + 1, 1] .- wfunc.points[1]
    z = wfunc.points[1, 1, cell_boundaries[3] + 1] .- wfunc.points[1]
    bloch = WannierFunction(wfunc.points, zeros(eltype(wfunc.values), size(wfunc.values)))
    dims = size(wfunc.values)
    for i1 in -3:1:3, i2 in -3:1:3, i3 in -3:1:3
        R_cryst = Vec3(i1, i2, i3)
        o1, o2, o3 = cell_boundaries .* R_cryst
        shiftvec = x * R_cryst[1] .+ y * R_cryst[2] .+ z * R_cryst[3]
        phase = ℯ^(2im * π * (R_cryst ⋅ kpoint))
        if i1 + i2 + i3 == 0
            continue
        end
        if i1 < 0
            o1 += i_neg_offset[1]
        elseif i1 > 0
            o1 += i_pos_offset[1]
        end
        if i2 < 0
            o2 += i_neg_offset[2]
        elseif i2 > 0
            o2 += i_pos_offset[2]
        end
        if i3 < 0
            o3 += i_neg_offset[3]
        elseif i3 > 0
            o3 += i_pos_offset[3]
        end

        for j3 in 1:dims[3]
            oid3 = mod1(j3 - o3, dims[3])
            for j2 in 1:dims[2]
                oid2 = mod1(j2 - o2, dims[2])
                for j1 in 1:dims[1]
                    oid1 = mod1(j1 - o1, dims[1])

                    bloch.values[j1, j2, j3] += phase * wfunc.values[oid1, oid2, oid3]
                end
            end
        end
    end
    return bloch
end

"Calculates the angular momentum between two wavefunctions and around the center."
function calc_angmom(
    wfc1::WannierFunction{N,T}, wfc2::WannierFunction{N,T}, center::Vec3{T}, cutoff=Inf
) where {N,T<:AbstractFloat}
    points = wfc1.points
    origin = points[1, 1, 1]
    da = points[2, 1, 1] - origin
    db = points[1, 2, 1] - origin
    dc = points[1, 1, 2] - origin
    V = SMatrix{3,3}(inv([convert(Array, da) convert(Array, db) convert(Array, dc)])')
    L = zero(Vec3{Complex{T}})
    c2 = cutoff^2
    @inbounds for i2 in 2:size(wfc1, 3)
        for i1 in 2:size(wfc1, 2)
            for i in 2:size(wfc1, 1)
                r = points[i, i1, i2] - center
                if dot(r, r) < c2
                    dw_cryst = Vec3(
                        wfc2.values[i, i1, i2] - wfc2.values[i - 1, i1, i2],
                        wfc2.values[i, i1, i2] - wfc2.values[i, i1 - 1, i2],
                        wfc2.values[i, i1, i2] - wfc2.values[i, i1, i2 - 1],
                    )

                    dw_cart = V * dw_cryst
                    L += (wfc1.values[i, i1, i2]',) .* cross(r, dw_cart)
                end
            end
        end
    end
    return -1im * L
end

#this doesn't work I think
# function calc_angmom_squared(wfc1::WannierFunction{N, T}, wfc2::WannierFunction{N, T}, center::Vec3{T}) where {N, T <: AbstractFloat}
#        points = wfc1.points
#     origin = points[1, 1, 1]
#     da     = points[2, 1, 1] - origin
#     db     = points[1, 2, 1] - origin
#     dc     = points[1, 1, 2] - origin
#     V      = SMatrix{3,3}(inv([convert(Array, da) convert(Array, db) convert(Array, dc)])')
#     Lsq      = zero(Complex{T})

#     @inbounds for i2 = 2:size(wfc1)[3]
#         for i1 = 2:size(wfc1)[2]
#             for i = 2:size(wfc1)[1]
#                 dw_cryst = Vec3(wfc2.values[i, i1, i2] - wfc2.values[i-1, i1,   i2],
#                                   wfc2.values[i, i1, i2] - wfc2.values[i,   i1-1, i2],
#                                                          wfc2.values[i, i1, i2] - wfc2.values[i,   i1,   i2-1])
#                    dw_cryst_sq = map(x->x .^2,dw_cryst)

#                 r       = points[i, i1, i2] - center
#                 dw_cart = V * dw_cryst
#                 Lsq    += wfc1.values[i, i1, i2] ⋅ (r[2]^2 * (dw_cryst_sq[1] + dw_cryst_sq[3]) +
#                                                     r[1]^2 * (dw_cryst_sq[2] + dw_cryst_sq[3]) +
#                                                     r[3]^2 * (dw_cryst_sq[1] + dw_cryst_sq[2]) -
#                                                     2 * (r ⋅ dw_cryst +
#                                                          r[2] * r[3] * dw_cryst[2] .* dw_cryst[3] +
#                                                          r[1] * r[3] * dw_cryst[1] .* dw_cryst[3] +
#                                                          r[1] * r[2] * dw_cryst[1] .* dw_cryst[2])) 

#                end
#            end
#        end
#        return Lsq
# end

function calc_spin(
    wfc1::WannierFunction{2,T}, wfc2::WannierFunction{2,T}
) where {T<:AbstractFloat}
    S = Vec3(
        SMatrix{2,2}(0, 1, 1, 0) / 2,
        SMatrix{2,2}(0, -1im, 1im, 0) / 2,
        SMatrix{2,2}(1, 0, 0, -1) / 2,
    )

    outS = zero(Vec3{Complex{T}})
    for (w1, w2) in zip(wfc1.values, wfc2.values)
        outS += (w1',) .* S .* (w2,)
    end
    return outS
end

"Calculates the dipole term between two wavefunctions. Make sure the wavefunctions are normalized!"
function calc_dip(
    wfc1::WannierFunction{N,T}, wfc2::WannierFunction{N,T}
) where {N,T<:AbstractFloat}
    out = zero(Vec3{Complex{T}})
    for (w1, w2, p) in zip(wfc1, wfc2, wfc1.points)
        out += w1' * w2 * p
    end
    return real(out)
end

# Is this code actually correct?
# "Calculates the dipoles from the supplied wannier dipole output."
# function calc_k_dips(dip_raw::Array{Tuple{Int, Int, Int, Int, Int, Vec3{T}}}, kpoints::AbstractArray) where T<:AbstractFloat
#        dim = 0
#        for i=1:length(dip_raw)
#               d = dip_raw[i][4]
#               if d > dim
#                      dim = d
#               else
#                      break
#               end
#        end
#        out = zeros(Vec3{T}, dim, dim)
#        tmp = [[zeros(Complex{T}, 3)] for i=1:dim, i1=1:dim]
#        for i=1:size(dip_raw)[1]
#               d = dip_raw[i]
#               complex_part = 2π*(kpoints[1]*d[1]+kpoints[2]*d[2]+kpoints[3]*d[3])
#               factor = exp(-1im * complex_part)
#               tmp[d[4],d[5]][1] += d[6][1] * factor
#               tmp[d[4],d[5]][2] += d[6][2] * factor
#               tmp[d[4],d[5]][3] += d[6][3] * factor
#        end
#        for i in eachindex(out)
#               out[i] = Vec3(real(tmp[i][1]),real(tmp[i][2]),real(tmp[i][3]))
#        end
#        return Mat{2*dim, 2*dim, Vec3{T}}([out zeros(out);zeros(out) out])
# end
