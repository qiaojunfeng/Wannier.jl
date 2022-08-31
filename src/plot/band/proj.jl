using ArgParse: ArgParse
import Wannier as Wan
import Plots as Pl

function parse_commandline(args=nothing)
    s = ArgParse.ArgParseSettings()
    ArgParse.@add_arg_table s begin
        # "--opt1"
        #     help = "an option with an argument"
        "--fermi_energy", "-f"
        help = "Fermi energy"
        arg_type = Float64
        default = 0.0
        "--show_orbitals"
        help = "Show each orbitals, e.g. s,p,d"
        action = :store_true
        "--dis_proj_min"
        help = "dis_proj_min"
        arg_type = Float64
        required = false
        "--dis_proj_max"
        help = "dis_proj_max"
        arg_type = Float64
        required = false
        "qebands"
        help = "Filename of QE bands.x output bands.dat file"
        required = true
        "qeprojs"
        help = "Filename of QE projwfc.x output prefix.proj.dat.projwfc_up file"
        required = true
    end
    if args == nothing
        return ArgParse.parse_args(s)
    else
        return ArgParse.parse_args(args, s)
    end
end

function main(args)
    parsed_args = parse_commandline(args)

    f_qe_bands = parsed_args["qebands"]
    f_qe_projs = parsed_args["qeprojs"]
    fermi_energy = parsed_args["fermi_energy"]
    show_orbitals = parsed_args["show_orbitals"]

    qe_bands = Wan.InOut.read_qe_bands(f_qe_bands)
    qe_projs = Wan.InOut.read_qe_projwfcup(f_qe_projs)

    # need to add labels
    qe_bands.symm_points_label = [
        "G", "X", "P", "N", "G", "M", "S", "S0", "G", "X", "R", "G", "M"
    ]
    #qe_bands.symm_points_label = ["L", "G", "X", "X", "G"]
    # fermi_energy = 1.5135500000e+01

    kwargs = Dict{Symbol,Any}()
    dis_proj_min = parsed_args["dis_proj_min"]
    dis_proj_max = parsed_args["dis_proj_max"]
    if dis_proj_min !== nothing && dis_proj_max !== nothing
        # only set the first atomic wfc
        colors = dropdims(sum(qe_projs.proj; dims=3); dims=3)
        qe_projs.proj .= 0.0
        p = view(qe_projs.proj, :, :, 1)
        p[colors .< dis_proj_min] .= 0.0
        p[colors .>= dis_proj_min] .= 0.5
        p[colors .>= dis_proj_max] .= 1.0
        kwargs[:color] = :jet1
        kwargs[:title] = "dis_proj_min/max = $(dis_proj_min)/$(dis_proj_max)"
        kwargs[:colorbar] = false
        kwargs[:marker_z] = p
    end

    #Pl.plotly()
    plt = Wan.plot_bands_projectabilities(
        qe_bands,
        qe_projs;
        fermi_energy=fermi_energy,
        show_orbitals=show_orbitals,
        show_gui=false,
        kwargs...,
    )
    # emin, emax = 15, 18
    # emin, emax = -8, 18
    # Pl.ylims!(plt, (emin, emax))
    Pl.gui(plt)

    # print("Hit <enter> to continue")
    # readline()
    print("Save figure to PDF? (Y/n)")
    y = lowercase(strip(readline()))
    if y == "y" || y == ""
        pdfname = "$(basename(f_qe_bands))+$(basename(f_qe_projs)).pdf"
        Pl.savefig(pdfname)
        println("Saved to $pdfname")
    end
end

main()

# using Colors

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
