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
            arg_type = Union{Int, Float64}
            default = 0
        # "--flag1"
        #     help = "an option without argument, i.e. a flag"
        #     action = :store_true
        "qe"
            help = "Filename of QE bands.x output bands.dat file"
            required = true
        "wannier90"
            help = "Seedname of Wannier90 win file"
            required = true
    end
    return ArgParse.parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    
    f_qe_bands = parsed_args["qe"]
    f_w90_bands = parsed_args["wannier90"]
    fermi_energy = parsed_args["fermi_energy"]

    qe_bands = Wan.InOut.read_qe_bands(f_qe_bands)
    w90_bands = Wan.InOut.read_wannier90_bands(f_w90_bands)

    # QE bands kpaths in `alat`, needs to rescale
    fac = (w90_bands.kpaths[end] - w90_bands.kpaths[1]) / (qe_bands.kpaths[end] - qe_bands.kpaths[1])
    qe_bands.kpaths .*= fac

    Wan.plot_bands_diff(qe_bands, w90_bands; fermi_energy=fermi_energy, label1="QE", label2="W90")

    # print("Hit <enter> to continue")
    # readline()
    print("Save figure to PDF? (Y/n)")
    y = lowercase(strip(readline()))
    if y == "y" || y == ""
        pdfname = "$(basename(f_qe_bands))+$(basename(f_w90_bands)).pdf"
        Pl.savefig(pdfname)
        println("Saved to $pdfname")
    end
end

main()
