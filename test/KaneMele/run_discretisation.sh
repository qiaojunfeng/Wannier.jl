#!/bin/sh

#Run the wannierize algorithm for different numbers of k-points.

system='kanemele'
method='parallelTransport'
# method='logInterpolation'
l_nu='0.0'
und='_'

for N in 25 50 75 100 #125 150 175 200
do
    python2 kanemele.py $N $l_nu
done


for N in 25 50 75 100 #125 150 175 200
do
    julia ../../run_wannierize.jl $system$und$l_nu$und$N $method >> results_$system$und$method.txt
done
