#!/bin/bash

set -e

NP=8
NK=4

PWX=pw.x
OGX=open_grid.x
W90X=wannier90.x
P2WX=pw2wannier90.x
BANDSX=bands.x

export OMP_NUM_THREADS=1

F='scf'
mpirun -n $NP $PWX -nk $NK -in $F.in > $F.out

mkdir bands
cp bands.in bandsx.in bands/
cd bands
cp -r ../out .
F='bands'
mpirun -n $NP $PWX -nk $NK -in $F.in > $F.out
F='bandsx'
mpirun -n $NP $BANDSX -nk $NK -in $F.in > $F.out
cd ..

# F='nscf'
# mpirun -n $NP $PWX -nk $NK -in $F.in > $F.out

F='opengrid'
# open_grid.x does not support kpool
mpirun -n $NP $OGX -in $F.in > $F.out

F='gaas.win'
mpirun -n $NP $W90X -pp $F

F='p2w'
mpirun -n $NP $P2WX -in $F.in > $F.out

F='gaas.win'
mpirun -n $NP $W90X $F
