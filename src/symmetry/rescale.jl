export rescale_repmat


function rescale_repmat(d)
    #not sure why all this is needed
    nbnd = size(d[1].h[1])[1]

    for ik in 1:size(d)[1]
        for isym in 1:size(d[ik].h)[1]
            for n in 1:nbnd
                vall = sum(abs.(d[ik].h[isym][n,:]))
                val  = abs(d[ik].h[isym][n,n])
                if (abs(vall-val) < 1.e-08) && !(round(val) == 0)
                    #if (abs(round(val)-val)>0.1) && (n<nbnd)
                    #end
                    d[ik].h[isym][n,n] = d[ik].h[isym][n,n]*round(val)/val 
                end 
            end
        end
    end
end 