#using Plots
#println("Creating animation of the obstruction")
#anim = @animate for i=1:2*p.N2-1
#               plot(real(Uint[:,abs(p.N2-i)+1,2,2]),imag(Uint[:,abs(p.N2-i)+1,2,2]),xlims=[-1,1],ylims=[-1,1])
#end

#gif(anim,"obs.gif",fps=25)
#close("all")
#Ucol1 = zeros(size(Uint)[1],size(Uint)[2],3)
#
#for i = 1:p.N1, j=1:p.N2
#	if Uint[i,j,2,1] != 0.0
#		phase = imag(log(im*Uint[i,j,2,1]/(abs(Uint[i,j,2,1]))))
#	else
#		phase = 0.0
#	end
#	array = Uint[i,j,:,1].*exp(-im*phase)
#
#	Ucol1[i,j,:] = [real(array[1]),imag(array[1]),imag(array[2])]
#end
