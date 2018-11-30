using PyCall
pygui(:tk)
using LinearAlgebra

include("wannierize.jl")
PyPlot.rc("font", family="serif")
PyPlot.rc("xtick", labelsize="x-small")
PyPlot.rc("ytick", labelsize="x-small")
PyPlot.rc("figure", figsize=(4,3))
#PyPlot.rc("text", usetex=false)


# Either run in batch mode, and take arguments from the CLI, or run in interactive mode and specify the file names manually
if(length(ARGS) >= 1)
    filename = ARGS[1]
    methodStr = ARGS[2]
    if methodStr == "parallelTransport"
	method = "parallel transport"
    elseif methodStr == "logInterpolation"
	method = "log interpolation"
    else
	println("ERROR: interpolation method not recognized.")
    end


    if length(ARGS) >= 3
	nbeg = parse(Int64,ARGS[5])
	nend = parse(Int64,ARGS[6])
    else
	nbeg = 1
	nend = Inf
    end
else
    filename = "tests/KaneMele/kanemele_0.0_25"
    #filename = "tests/silicon/silicon"
    method = "parallel transport"
    #method = "log interpolation"
end

read_amn = false #true
read_eig = false #true

p = read_system(filename,read_amn,read_eig)

if(length(ARGS) <1)
    nbeg = 1
    nend = p.nwannier
end

println("Method of interpolation: $method")

println("$(p.nband) bands, wannierizing bands $nbeg:$nend, $(p.N1) x $(p.N2) x $(p.N3) grid, $(p.nntot) neighbors")

#Build the gauge that makes the Bloch frame continuous on the Brillouin Zone.
#This is equivalent to building a set of algebraic decaying Wannier functions
p, interp = make_wannier(p,method)

#Plot and display the results
close("all")
plot_results(p,interp)
print_error(p)
