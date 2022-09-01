#!/bin/bash

set -e

NP=8
NK=8

PWX=pw.x
W90X=wannier90.x
P2WX=pw2wannier90.x

F='scf'
mpirun -n $NP $PWX -nk $NK -in $F.in > $F.out

F='nscf'
mpirun -n $NP $PWX -nk $NK -in $F.in > $F.out

F='graphene.win'
mpirun -n $NP $W90X -pp $F

F='p2w'
mpirun -n $NP $P2WX -in $F.in > $F.out

F='graphene.win'
mpirun -n $NP $W90X $F
