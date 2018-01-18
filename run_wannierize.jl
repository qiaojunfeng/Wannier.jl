include("wannierize.jl")

filename = "silicon"

read_amn = false
read_eig = false

p = read_system(filename,read_amn,read_eig)

nbeg = 1
nend = p.nwannier

#Set the seed of the pseudo-random generator
#srand(0)

println("$(p.nband) bands, wannierizing bands $nbeg:$nend, $(p.N1) x $(p.N2) x $(p.N3) grid, $(p.nntot) neighbors")

#Build the gauge that makes the Bloch frame continuous on the Brillouin Zone.
#This is equivalent to building a set of algebraic decaying Wannier functions
p, Obs_array, Uint = make_wannier(p)

# #Plot and display the results
# plot_results(p)
# print_error(p)

# figure()
# for i=1:5:p.N2
#     plot(real(Uint[:,i,2,2]),imag(Uint[:,i,2,2]),"x-")
# end
# title("real 22 / imag 22")

# figure()
# for i=1:10:p.N2
#     plot(real(Uint[:,i,1,2]),imag(Uint[:,i,1,2]),"x-")
# end
# title("real 12 / imag 12")

# figure()
# for i=1:5:p.N2
#     plot(real(Uint[:,i,:,2]))
# end

# #using Plots
# #println("Creating animation of the obstruction")
# #anim = @animate for i=1:2*p.N2-1
# #               plot(real(Uint[:,abs(p.N2-i)+1,2,2]),imag(Uint[:,abs(p.N2-i)+1,2,2]),xlims=[-1,1],ylims=[-1,1])
# #end
# #
# #gif(anim,"obs.gif",fps=25)
