using LinearAlgebra
using Brillouin

"""
Get kpoint coordinates from a KPath, same as W90.

Use the kpath density of first segment to generate the following kpaths,
also need to take care of high symmetry kpoints at the start and end of each segment.

n_points: number of kpoints in the first segment, remaining segments
    have the same density as 1st segment.
"""
function interpolate_w90(kpath::KPath, n_points::Int)
    # cartesian
    kpath_cart = cartesianize(kpath)
    # kpath density from first two kpoints
    k1, k2 = kpath_cart.paths[1][1:2]
    seg = kpath_cart.points[k2] - kpath_cart.points[k1]
    seg_norm = norm(seg)
    dk = seg_norm / n_points

    # kpoints along path
    kpaths = Vector{Vector{Vec3{Float64}}}()
    # symmetry points
    labels = Vector{Dict{Int,Symbol}}()

    for path in kpath_cart.paths
        kpaths_line = Vector{Vec3{Float64}}()
        labels_line = Dict{Int,Symbol}()

        n_seg = length(path) - 1
        n_x_line = 0
        for j in 1:n_seg
            k1 = path[j]
            k2 = path[j + 1]

            seg = kpath_cart.points[k2] - kpath_cart.points[k1]
            seg_norm = norm(seg)

            n_x_seg = Int(round(seg_norm / dk))
            x_seg = collect(range(0, seg_norm, n_x_seg + 1))
            dvec = seg / seg_norm

            # column vector * row vector = matrix
            kpt_seg = dvec * x_seg'
            kpt_seg .+= kpath_cart.points[k1]

            if j == 1
                push!(labels_line, 1 => k1)
            else
                # remove repeated points
                popfirst!(x_seg)
                kpt_seg = kpt_seg[:, 2:end]
            end
            n_x_line += length(x_seg)
            push!(labels_line, n_x_line => k2)

            append!(kpaths_line, [v for v in eachcol(kpt_seg)])
        end

        push!(kpaths, kpaths_line)
        push!(labels, labels_line)
    end

    basis = kpath.basis
    setting = Ref(Brillouin.CARTESIAN)
    kpi = KPathInterpolant(kpaths, labels, basis, setting)
    # to fractional
    latticize!(kpi)

    return kpi
end

"""
Get x axis value for plotting, cartesian length
"""
function get_x(kpi::KPathInterpolant)
    kpi_cart = cartesianize(kpi)
    x = Vector{Float64}()

    for path in kpi_cart.kpaths
        n_points = length(path)

        push!(x, 0)
        for j in 2:n_points
            k1 = path[j - 1]
            k2 = path[j]
            dx = norm(k2 - k1)
            push!(x, dx)
        end
    end

    return cumsum(x)
end
