Julia code for the computation of Wannier functions.

# Functionalities implemented
- [1] E. Cancès, A. Levitt, G. Panati, G. Stoltz, "Robust determination of maximally-localized Wannier functions", https://arxiv.org/abs/1605.07201

(get Wannier functions for isolated bands without specifying an initial gauge)

- [2] A. Damle, A. Levitt, L. Lin, "Variational formulation for Wannier functions with entangled band structure", link to come

(alternative minimization algorithm to Wannier90 for both isolated and non-isolated bands)

# Usage
This assumes you are familiar with the Wannier90 workflow.

For the method in [1], write the `.amn/win` and use `run_wannierize.jl`. See `tests/silicon/run_all.sh` for the Wannierization of the first four bands of Silicon (where using a random initial guess fails to compute good Wannier functions but the method in [1] does).

For the method in [2], write the `.amn/mmn/eig/win` files and use `run_optim.jl` (see parameters in that file). See `tests/free` for the computation of Wannier functions for the free electron gas in [2].

Requirements: Julia (tested on v0.6), the Optim library for the method in [2], PyPlot for the plots, and wannier90/Quantum Espresso/Python for the tests.

# Contact
This is research code, not necessarily user-friendly or extremely robust. If you have questions or encounter problems, get in touch!
