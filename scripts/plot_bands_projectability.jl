#!/usr/bin/env julia
import ArgParse
import Wannier as Wan
import Plots as Pl

function parse_commandline()
    s = ArgParse.ArgParseSettings()
    @ArgParse.add_arg_table s begin
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
    return ArgParse.parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    
    f_qe_bands = parsed_args["qebands"]
    f_qe_projs = parsed_args["qeprojs"]
    fermi_energy = parsed_args["fermi_energy"]
    show_orbitals = parsed_args["show_orbitals"]

    qe_bands = Wan.InOut.read_qe_bands(f_qe_bands)
    qe_projs = Wan.InOut.read_qe_projwfcup(f_qe_projs)

    # need to add labels
    # qe_bands.symm_points_label = ["G", "X", "P", "N", "G", "M", "S", "S0", "G", "X", "R", "G", "M"]
    qe_bands.symm_points_label = ["L", "G", "X", "X", "G"]
    # fermi_energy = 1.5135500000e+01

    kwargs = Dict{Symbol,Any}()
    dis_proj_min = parsed_args["dis_proj_min"]
    dis_proj_max = parsed_args["dis_proj_max"]
    if dis_proj_min !== nothing && dis_proj_max !== nothing
        # only set the first atomic wfc
        colors = dropdims(sum(qe_projs.proj, dims=3), dims=3)
        qe_projs.proj .= 0.0
        p = view(qe_projs.proj,:,:,1)
        p[colors .< dis_proj_min] .= 0.0
        p[colors .>= dis_proj_min] .= 0.5
        p[colors .>= dis_proj_max] .= 1.0
        kwargs[:color] = :jet1
        kwargs[:title] = "dis_proj_min/max = $(dis_proj_min)/$(dis_proj_max)"
        kwargs[:colorbar] = false
    end

    #Pl.plotly()
    plt = Wan.plot_bands_projectabilities(qe_bands, qe_projs; 
        fermi_energy=fermi_energy, show_orbitals=show_orbitals, show_gui=false, kwargs...)
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
