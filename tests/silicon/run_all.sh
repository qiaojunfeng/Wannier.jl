#!/bin/sh

system=silicon
method=parallelTransport
#method=logInterpolation
#N=$1
N=5

#mkdir $method.$N/
cd $method.$N/


cp ../$system.win.template $system.win
echo mp_grid $N $N $N >> $system.win
echo "begin kpoints" >> $system.win
python2 ../change_kpts.py $N >> $system.win
echo "end kpoints" >> $system.win

cp ../$system.scf.template $system.scf
echo $N $N $N 0 0 0 >> $system.scf

cp ../$system.nscf.template $system.nscf
echo "$(($N * $N * $N))" >> $system.nscf
python2 ../change_kpts.py $N print_weights >> $system.nscf

echo "SCF"
pw.x < $system.scf > scf.out
echo "Non-SCF"
pw.x < $system.nscf > nscf.out
echo "wannier -pp"
wannier90.x -pp $system
echo "pw2wannier mmn"
pw2wannier90.x < ../$system.pw2wan.mmn > pw2wan.mmn.out
echo "pw2wannier amn"
pw2wannier90.x < ../$system.pw2wan.amn > pw2wan.amn.out
##echo "wannier90: plot bands before wannierisation"
##wannier90.x $system
##mkdir bands_before/
##mv $system\_band.* bands_before/
#
echo "wannierize"
julia ../../../run_wannierize.jl $system $method >> $system.wannierize.out

#TO PLOT BANDS
#cp ../$system.plotbands.win.template $system.win
#echo mp_grid $N $N $N >> $system.win
#echo "begin kpoints" >> $system.win
#python2 ../change_kpts.py $N >> $system.win
#echo "end kpoints" >> $system.win

echo "wannier90: compute spreads"
wannier90.x $system
#mkdir bands_after/
#mv $system\_band.* bands_after/
mkdir spreads/
grep SPRD $system.wout | cut -c 60-70 > spreads/$system.sprd
cd ../
