include("wannierize_utils.jl")

using PyPlot

## Assumptions: the kpoints are contained in a NxNxN cartesian grid, the neighbor list must contain the six cartesian neighbors

## N1xN2xN3 grid, filename.mmn must contain the overlaps
## nbeg and nend specify the window to wannierize
## Input Mmn file is nband x nband x nkpt x nntot
## Output is (nwannier = nend - nbeg + 1) x nband x nkpt, padded with zeros
function make_wannier(p)
    t1 = collect(0:p.N1-1)/p.N1
    t2 = collect(0:p.N2-1)/p.N2
    t3 = collect(0:p.N3-1)/p.N3
    Ntot = p.N1*p.N2*p.N3
    #nwannier = nend-nbeg+1

    A0 = deepcopy(p.A)
    M0 = deepcopy(p.M)
    A = p.A
    M = p.M
    neighbors = p.neighbors

    fill!(A,NaN) #protection: A must be filled by the algorithm
    println("Filling (k,0,0)")
    A[:,:,:,1,1] = propagate(eye(p.nwannier), [p.ijk_to_K[i,1,1] for i=1:p.N1], p)
    
                                                                                            
    # compute obstruction matrix                                                            
    Obs = normalize_matrix(overlap_A([p.N1,1,1],[1,1,1],p))
    #println("Obstruction matrix = ")                                                        
    #println(Obs)                                                                            
    d,V = eig(Obs)                                                                          
                                                                                            
    # and pull it back                                                                      
    for i=1:p.N1                                                                              
	A[:,:,i,1,1] = A[:,:,i,1,1]*powm(Obs,t1[i])                                         
    end
    
    println("Filling (k1,k2,0)")
    for i=1:p.N1
	A[:,:,i,:,1] = propagate(A[:,:,i,1,1], [p.ijk_to_K[i,j,1] for j=1:p.N2],p)
    end
    
    
    # corner obstruction
    Obs = normalize_matrix(overlap_A([1,p.N2,1],[1,1,1],p))
    d,V = eig(Obs)
    #println("Corner obstruction")
    #println(Obs)

    logd = log.(d)
    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i=1:p.nwannier
        if imag(logd[i]) < -pi+.01
            logd[i] = logd[i] + 2pi*im
        end
    end
    # pull it back
    for i=1:p.N1
        for j=1:p.N2
            A[:,:,i,j,1] = A[:,:,i,j,1]*V*diagm(exp.(t2[j]*logd))*V'
        end
    end
    
    # Pull back the line obstruction
    phases = zeros(p.N1,p.nwannier)
    Obs_array = zeros(Complex128,p.N1,p.nwannier,p.nwannier) 
    detObs = zeros(Complex128,p.N1)
    eigs = zeros(Complex128,p.N1,p.nwannier)
    for i =1:p.N1
    Obs_array[i,:,:] =  normalize_matrix(overlap_A([i,p.N2,1],[i,1,1],p)) 
    	detObs[i] = det(Obs_array[i,:,:])
    end
    
    # Find a continuous log of the determinant
    logDet = imag(log.(detObs))
    for i=2:p.N1
    	kmin = indmin([abs(logDet[i]+2*pi*k-logDet[i-1]) for k in -1:1])
    	logDet[i] = logDet[i]+(kmin-2)*2*pi
    end
    for i=1:p.N1
	#Obs_array[i,:,:] = [detObs[i]' 0; 0 1] * Obs_array[i,:,:]
	Obs_array[i,:,:] = exp(-im*logDet[i]/2) * Obs_array[i,:,:]
	eigs[i,:] = eigvals(Obs_array[i,:,:])
    end

    # Interpolate the line obstruction
    logd = Complex128[0; 0]
    Uint = zeros(Complex128,p.N1,p.N2,p.nwannier,p.nwannier)
    #println("matrix_pole = $matrix_pole")

    for i=1:p.N1
	Obs = Obs_array[i,:,:]
	#println("Obs_array[$i,:,:] = $(Obs_array[i,:,:])")
	d,V = eig(Obs)
	for j=1:p.N2
	    Uint[i,j,:,:] = powm(Obs,t2[j])
	end
	
        for j=1:p.N2
		#A[:,:,i,j,1] = A[:,:,i,j,1]*[exp(im*logDet[i]*t2[j]) 0; 0 1]
		A[:,:,i,j,1] = A[:,:,i,j,1] * exp(im*logDet[i]*t2[j]/2)
		A[:,:,i,j,1] = A[:,:,i,j,1] * Uint[i,j,:,:]
        end
    end

    # Propagate along the third dimension   
    println("Filling (k1,k2,k3)")
    for i=1:p.N1,j=1:p.N2
    	A[:,:,i,j,:] = propagate(A[:,:,i,j,1], [p.ijk_to_K[i,j,k] for k=1:p.N3],p)
    end
    
    # Fix corner
    Obs = normalize_matrix(overlap_A([1,1,p.N3],[1,1,1],p))
    d,V = eig(Obs)
    logd = imag(log.(d))
    #println("Obstruction matrix  - Id = $(norm(Obs-I)))")

    # Hack to avoid separating eigenvalues at -1. TODO understand that
    for i =1:p.nwannier
        if imag(logd[i]) < -pi+.01
            logd[i] = logd[i] + 2pi*im
        end
    end

    for k=1:p.N3
        # fixer = powm(Obs,t3[k])
        fixer = V*diagm(exp.(t3[k]*logd))*V'
        for i=1:p.N1,j=1:p.N2
            A[:,:,i,j,k] = A[:,:,i,j,k]*fixer
        end
    end
    
    # Fix first edge
    for i=1:p.N1
	    Obs = normalize_matrix(overlap_A([i,1,p.N3], [i,1,1],p))
        for k=1:p.N3
            fixer = powm(Obs, t3[k])
            for j=1:p.N3
                A[:,:,i,j,k] = A[:,:,i,j,k]*fixer
            end
        end
    end
    
    # Fix second edge
    for j=1:p.N2
	    Obs = normalize_matrix(overlap_A([1,j,p.N3], [1,j,1], p))
        for k=1:p.N3
            fixer = powm(Obs, t3[k])
            for i=1:p.N1
                A[:,:,i,j,k] = A[:,:,i,j,k]*fixer
            end
        end
    end

    # Fix whole surface
    for i=1:p.N1,j=1:p.N2
	    Obs = normalize_matrix(overlap_A([i,j,p.N3],[i,j,1],p)) 
        for k=1:p.N3
            A[:,:,i,j,k] = A[:,:,i,j,k]*powm(Obs,t3[k])
        end
    end

    #plot_Bloch_frame_slice(p,A0)

    write_amn(p,A,filename)

    return(p,Obs_array,Uint)
end

function plot_results(p)
    suffix = ""

    #Initialize arrays
    omegaright = zeros(p.N1,p.N2)
    omegaup = zeros(p.N1,p.N2)
    Aright1 = zeros(Complex128,p.N1,p.N2)
    Aright2 = zeros(Complex128,p.N1,p.N2)
    detA = zeros(Complex128,p.N1,p.N2)

    #Fill arrays
    imax = 1
    jmax = 1
    omegamax = 0.
    for i=1:p.N1,j=1:p.N2
        right = i==p.N1 ? 1 : i+1
        up = j==p.N2 ? 1 : j+1
        Aright1[i,j] = (normalize_matrix(overlap_A([right,j,1],[i,j,1],p)))[1,1]
        Aright2[i,j] = (normalize_matrix(overlap_A([right,j,1],[i,j,1],p)))[2,2]
        detA[i,j] = det(normalize_matrix(overlap_A([right,j,1],[i,j,1],p)))
        omegaright[i,j] = norm(normalize_matrix(overlap_A([right,j,1],[i,j,1],p)) - I)
        omegaup[i,j] = norm(normalize_matrix(overlap_A([i,up,1],[i,j,1],p)) - I)
        if (omegaright[i,j] + omegaup[i,j] > omegamax)
       	 omegamax = omegaright[i,j] + omegaup[i,j]
       	 imax = i
       	 jmax = j
        end

    end
    println("omegamax = $omegamax, [i,j] = $([imax,jmax])")

    matshow(omegaright)
    colorbar()
    title("Omega right N1,N2,N3 = $(p.N1), $(p.N2), $(p.N3)  $(p.filename) $suffix")
    savefig("omega_right_$(p.filename).$suffix.pdf")
    matshow(omegaup)
    colorbar()
    title("Omega up N1,N2,N3 = $(p.N1), $(p.N2), $(p.N3) $(p.filename) $suffix")
    savefig("omega_up_$(p.filename).$suffix.pdf")

    plot_surface_obstructions(p, "_1_none")
    plot_surface_obstructions(p, "_2_corners")
    plot_surface_obstructions(p, "_3_edges")
    plot_surface_obstructions(p, "_4_surface")
end

function print_error(p)
    err_before = 0.
    err_after  = 0.
    for i=1:p.N1, j=1:p.N2, k=1:p.N3
	    err_before += norm(normalize_matrix(overlap(p.ijk_to_K[i,j,k],p.ijk_to_K[(i%p.N1+1),j,k],p)) - I)^2
	    err_after  += norm(normalize_matrix(overlap_A([i,j,k],[(i%p.N1+1),j,k],p))- I)^2
	    err_before += norm(normalize_matrix(overlap(p.ijk_to_K[i,j,k],p.ijk_to_K[i,(j%p.N2)+1,k],p)) - I)^2
	    err_after  += norm(normalize_matrix(overlap_A([i,j,k],[i,(j%p.N2)+1,k],p)) - I)^2
    end
    err_before =  sqrt( 1/(p.N1*p.N2*p.N3) * err_before)
    err_after =  sqrt(1/(p.N1*p.N2*p.N3) * err_after)
    println("err before = $(err_before)")
    println("err after  = $(err_after)") 
end

function plot_Bloch_frame_slice(p,A0)
    if p.N3 > 2
    	N3_list = [1, div(p.N3,2), p.N3]
    else
	N3_list = [1]
    end

    for slice in N3_list
	matshow(real(A0[slice,1,1,:,:]))
	colorbar()
	title("Bloch frame slice before algorithm, N3 = $slice")
	savefig("bloch_frame_before_$slice.pdf")
    end
    for slice in N3_list
	matshow(real(p.A[slice,1,1,:,:]))
	colorbar()
	title("Bloch frame slice after algorithm, N3 = $slice")
	savefig("bloch_frame_after_$slice.pdf")
    end


end



