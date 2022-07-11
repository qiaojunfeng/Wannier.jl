import Plots as Pl
using Colors: Colors

function merge_consecutive_labels(
    symm_idx::AbstractArray{Int}, symm_label::AbstractVector{String}
)
    idx = copy(symm_idx)
    labels = copy(symm_label)

    counter = 2
    for i in 2:length(symm_idx)
        if symm_idx[i] == symm_idx[i - 1] + 1
            labels[counter - 1] *= "|$(symm_label[i])"
        else
            idx[counter] = symm_idx[i]
            labels[counter] = symm_label[i]
            counter += 1
        end
    end

    return idx[1:(counter - 1)], labels[1:(counter - 1)]
end

function plot_band(x::AbstractVector{T}, E::AbstractArray{T}; kwargs...) where {T<:Real}
    Pl.PlotlyBackend()

    plt = Pl.plot()

    plot_band!(plt, x, E; kwargs...)

    return plt
end

function plot_band!(
    plt,
    x::AbstractVector{T},
    E::AbstractArray{T};
    fermi_energy::Union{Nothing,T}=nothing,
    symm_idx::Union{Nothing,AbstractArray{Int}}=nothing,
    symm_label::Union{Nothing,AbstractVector{String}}=nothing,
    kwargs...,
) where {T<:Real}
    # color=, line_z::Union{Matrix{Float64},Missing}=missing)

    ndims(E) > 2 && error("E must be a 1D or 2D array")

    Pl.ylabel!(plt, "Energy (eV)")

    if fermi_energy !== nothing
        Pl.hline!(plt, [fermi_energy]; linestyle=:dash, linecolor=:blue, linewidth=0.2)
    end

    if symm_idx !== nothing
        n_symm = length(symm_idx)
        n_symm != length(symm_label) &&
            error("symm_idx and symm_label must have the same length")

        idx, label = merge_consecutive_labels(symm_idx, symm_label)
        Pl.vline!(plt, x[idx]; linecolor=:black, linewidth=0.2)
        # decrease a bit, otherwise the last label is not shown
        xtick = x[idx]
        xtick[end] -= 0.5 * (x[end] - x[end - 1])
        Pl.xticks!(plt, xtick, label)
    end

    # # kwargs is immutable, copy to a new Dict
    # varargs = Dict{Symbol,Any}()
    # for (k, v) in kwargs
    #     varargs[k] = v
    # end
    # @debug "varargs" varargs
    # if !haskey(varargs, :color)
    #     varargs[:color] = :copper
    # end

    E_plot = E
    if ndims(E) == 2
        E_plot = E'
    end

    Pl.plot!(plt, x, E_plot; legend=false, grid=false, framestyle=:box, kwargs...)

    xmin = x[1]
    xmax = x[end]
    Pl.xlims!(plt, (xmin, xmax))

    emin = minimum(E) - 0.5
    emax = maximum(E) + 0.5
    Pl.ylims!(plt, (emin, emax))

    return nothing
end

# function plot_band_projectabilities(bands::Bands, projectabilities::Projectabilities;
#     fermi_energy::Union{Int,Float64}=0.0, show_orbitals::Bool=false, show_gui=true, kwargs...)

#     plt = nothing
#     if !show_orbitals
#         colors = dropdims(sum(projectabilities.proj, dims=3), dims=3)
#         plt = plot_band(Pl.plot(), bands; fermi_energy=fermi_energy, line_z=colors, kwargs...)
#     else
#         # merge atomic wfcs, e.g. px+py+pz -> p
#         wfc_typs = Dict()
#         for iw = 1:projectabilities.num_wfcs
#             lab = projectabilities.wfcs_type[iw].atom_label * projectabilities.wfcs_type[iw].wfc_label
#             if haskey(wfc_typs, lab)
#                 wfc_typs[lab] += projectabilities.proj[:, :, iw]
#             else
#                 wfc_typs[lab] = projectabilities.proj[:, :, iw]
#             end
#         end
#         num_subplots = length(wfc_typs)
#         @debug "wfc_typs" wfc_typs
#         # dis_colors = Colors.distinguishable_colors(num_subplots)
#         layout = Pl.grid(1, num_subplots)
#         plots = []
#         clims = (0.0, 1.0)
#         for typ in keys(wfc_typs)
#             # p = Pl.plot(bands.kpaths, bands.energies; line_z=colors, color=dis_colors[iw],
#             # legend=false, grid=false, framestyle=:box)
#             p = plot_band(Pl.plot(), bands;
#                 fermi_energy=fermi_energy, line_z=wfc_typs[typ],
#                 # color=dis_colors[iw],
#                 # clims=clims,
#                 colorbar=true,
#                 kwargs...)
#             Pl.title!(p, typ)
#             push!(plots, p)
#         end
#         plt = Pl.plot(plots..., layout=layout, link=:all)
#     end

#     if show_gui
#         Pl.gui(plt)
#     end

#     return plt
# end

# function plot_bands_diff(bands1::Bands, bands2::Bands; fermi_energy::Union{Int,Float64}=0.0,
#     dis_froz_max::Union{Nothing,Float64}=nothing, label1::String="Bands1", label2::String="Bands2")
#     plt = Pl.plot(grid=false, framestyle=:box)

#     Pl.ylabel!(plt, "Energy (eV)")
#     Pl.hline!(plt, [fermi_energy]; linestyle=:dash, linecolor=:blue, linewidth=0.2, label="Fermi energy")

#     if dis_froz_max !== nothing
#         Pl.hline!(plt, [dis_froz_max]; linestyle=:dash, linecolor=:green, linewidth=0.2, label="dis_froz_max")
#     end

#     if bands1.num_symm_points != 0 && all(bands1.symm_points_label .!= "")
#         merged_labels = merge_consecutive_labels(bands1.symm_points, bands1.symm_points_label)
#         Pl.vline!(plt, bands1.kpaths[bands1.symm_points]; linecolor=:black, linewidth=0.2, label="")
#         Pl.xticks!(plt, bands1.kpaths[bands1.symm_points], merged_labels)
#     end
#     if bands2.num_symm_points != 0 && all(bands2.symm_points_label .!= "")
#         merged_labels = merge_consecutive_labels(bands2.symm_points, bands2.symm_points_label)
#         Pl.vline!(plt, bands2.kpaths[bands2.symm_points]; linecolor=:black, linewidth=0.2, label="")
#         Pl.xticks!(plt, bands2.kpaths[bands2.symm_points], merged_labels)
#     end

#     # only show 1 label in the legend
#     labels = hcat(label1, fill("", 1, bands1.num_bands))
#     Pl.plot!(plt, bands1.kpaths, bands1.energies; linecolor=:black, label=labels)
#     labels = hcat(label2, fill("", 1, bands2.num_bands))
#     Pl.plot!(plt, bands2.kpaths, bands2.energies; linecolor=:red, linestyle=:dash, label=labels)

#     xmin = min(bands1.kpaths[1], bands2.kpaths[1])
#     xmax = max(bands1.kpaths[end], bands2.kpaths[end])
#     Pl.xlims!(plt, (xmin, xmax))

#     emin = min(minimum(bands1.energies), minimum(bands2.energies)) - 0.5
#     emax = max(maximum(bands1.energies), maximum(bands2.energies)) + 0.5
#     Pl.ylims!(plt, (emin, emax))

#     Pl.gui(plt)

#     return plt
# end
