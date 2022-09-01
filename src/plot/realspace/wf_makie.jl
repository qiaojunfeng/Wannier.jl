using LinearAlgebra
using GeometryBasics
using GLMakie

# function plot_wf_makie(
#     W::AbstractArray{T,3};
#     iso::Union{T,Nothing}=nothing,
# ) where {T<:Real,U<:Integer}
#     scene = Makie.Scene()
#     cam3d!(scene)

#     if isnothing(iso)
#         iso = _guess_isolevel(W)
#     end
#     a, b = minimum(W), maximum(W)

#     if iso >= a
#         algo = MarchingCubes(; iso=iso, insidepositive=true)
#         mc = GeometryBasics.Mesh(W, algo)
#         Makie.mesh!(scene, mc; color="#C3423F", shading=true)#color=[norm(v) for v in coordinates(mc)])
#     end
#     if iso <= b
#         algo = MarchingCubes(; iso=-iso, insidepositive=true)
#         mc = GeometryBasics.Mesh(W, algo)
#         Makie.mesh!(scene, mc; color="#5BC0EB", shading=true)#color=[norm(v) for v in coordinates(mc)])
#     end

#     return scene
# end
