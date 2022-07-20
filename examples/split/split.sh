#!/bin/bash
# Wrapper to call w90chk2chk.x, split Wannierization, and prepare win files.

# exit immediately if a command exits with a non-zero status.
set -e

# change these to your own paths
w90chk2chk="$HOME/git/wannier90/w90chk2chk.x"
wannier="$HOME/.julia/bin/wannier"

valence_dir='val'
conduction_dir='cond'

# pass at least a seedname
if [ $# -lt 1 ]; then
    echo "Usage: $0 <seedname>"
    exit 1
fi
# the last argument is the seedname
seedname="${@: -1}"

# From binary to textual
$w90chk2chk -export $seedname
printf '\n\n'

# split Wannierization
$wannier splitvc $@
printf '\n\n'

# function to prepare win file
modify_win() {
    local seedname=$1
    # modify num_bands, num_wann
    local line=$(head -n 2 $seedname.amn | tail -n 1)
    local num_bands=$(echo $line | cut -d' ' -f1)
    local num_wann=$(echo $line | cut -d' ' -f3)
    sed -i "s/num_bands[[:space:]]*=.*/num_bands = $num_bands/" $seedname.win
    sed -i "s/num_wann[[:space:]]*=.*/num_wann = $num_wann/" $seedname.win
}

# Valence
cd "$valence_dir"

if test ! -e "$seedname.win"; then
    cp ../$seedname.win .
    modify_win $seedname
fi

cd ../
printf '\n\n'

# Conduction
cd "$conduction_dir"

if test ! -e "$seedname.win"; then
    cp ../$seedname.win .
    modify_win $seedname
fi

cd ../
printf '\n\n'
