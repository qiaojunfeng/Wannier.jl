### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 8d696b84-fc3f-11ec-340f-c3538d107d41
using Pkg

# ╔═╡ d07812ac-fbde-4f20-8dba-21d89e689f12
Pkg.activate("/home/junfeng/git/Wannier.jl")

# ╔═╡ c4f88233-e529-4ca2-9a86-a97a0bb072b2
using Revise

# ╔═╡ 065e0a97-d57d-4322-a09d-befa7a0204fe
using Plots

# ╔═╡ 06961146-7fa6-42a0-b922-423cbe3a64c0
using Wannier

# ╔═╡ 9ae5af97-7c12-452a-9075-b208b867490e
plotly();

# ╔═╡ 64cada09-26fa-437b-bd54-fe733e12aa8d
const NOTEBOOK_PATH = replace(@__FILE__, r"#==#.*" => "")

# ╔═╡ c383427b-9ad4-4047-adfe-b65cbea5d9aa
const DATA_PATH = joinpath(
    splitpath(NOTEBOOK_PATH)[1:(end - 3)]..., "test/fixtures/valence/band"
)

# ╔═╡ 59835aaa-7b53-4b50-8cdf-0bc21a504a3b
cd(DATA_PATH)

# ╔═╡ 2ef8c7d5-a64a-40c6-a673-9368108c1cbb
readdir()

# ╔═╡ 39c610b0-3634-4341-8eb2-aa8a4722dc0c
band = read_w90_bands("silicon")

# ╔═╡ 49ebb057-d25b-45e8-92fa-9b074a80aa94
Wannier.plot_band(band.x, band.E; symm_idx=band.symm_idx, symm_label=band.symm_label)

# ╔═╡ 8315073f-f96f-4fa7-9ae4-362c38b73ac0
model = read_seedname("silicon");

# ╔═╡ 770646d5-daa0-4afb-b37e-c9b98e942824
E = Wannier.interpolate(model, band.kpoints)

# ╔═╡ 13c7cef6-4f92-454e-9504-2dfe70753a9a
Wannier.plot_band(band.x, E; symm_idx=band.symm_idx, symm_label=band.symm_label)

# ╔═╡ Cell order:
# ╠═8d696b84-fc3f-11ec-340f-c3538d107d41
# ╠═d07812ac-fbde-4f20-8dba-21d89e689f12
# ╠═c4f88233-e529-4ca2-9a86-a97a0bb072b2
# ╠═065e0a97-d57d-4322-a09d-befa7a0204fe
# ╠═9ae5af97-7c12-452a-9075-b208b867490e
# ╠═06961146-7fa6-42a0-b922-423cbe3a64c0
# ╠═64cada09-26fa-437b-bd54-fe733e12aa8d
# ╠═c383427b-9ad4-4047-adfe-b65cbea5d9aa
# ╠═59835aaa-7b53-4b50-8cdf-0bc21a504a3b
# ╠═2ef8c7d5-a64a-40c6-a673-9368108c1cbb
# ╠═39c610b0-3634-4341-8eb2-aa8a4722dc0c
# ╠═49ebb057-d25b-45e8-92fa-9b074a80aa94
# ╠═8315073f-f96f-4fa7-9ae4-362c38b73ac0
# ╠═770646d5-daa0-4afb-b37e-c9b98e942824
# ╠═13c7cef6-4f92-454e-9504-2dfe70753a9a
