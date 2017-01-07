Companion code for https://arxiv.org/abs/1605.07201

Needs julia (v0.5) and PyPlot for the plots, and quantum espresso/python for the tests. Usage: compute system.mmn with pw2wannier and then run julia wannierize.jl system. The k-point mesh is assumed to be cartesian, with the third dimension as fast index (see tests).
