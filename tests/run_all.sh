#!/bin/sh

system=silicon
N=10
WANNIER90=~/wannier90/wannier90.x

cp $system.win.template $system.win
echo mp_grid $N $N $N >> $system.win
echo "begin kpoints" >> $system.win
python change_kpts.py $N >> $system.win
echo "end kpoints" >> $system.win

cp $system.nscf.template $system.nscf
echo "$(($N * $N * $N))" >> $system.nscf
python change_kpts.py $N print_weights >> $system.nscf

echo "SCF"
pw.x < $system.scf > scf.out
echo "Non-SCF"
pw.x < $system.nscf > nscf.out
echo "wannier -pp"
$WANNIER90 -pp $system
echo "pw2wannier mmn"
pw2wannier90.x < $system.pw2wan.mmn > pw2wan.mmn.out
echo "pw2wannier amn"
pw2wannier90.x < $system.pw2wan.amn > pw2wan.amn.out
echo "wannierize"
julia wannierize.jl
echo "wannier90"
$WANNIER90 $system
grep SPRD $system.wout | cut -c 60-70 > spreads/$N/$system.sprd
