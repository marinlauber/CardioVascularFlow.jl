using StaticArrays
using WaterLily
using MeshIO
using StaticArrays
using WriteVTK
using NearestNeighbors

include("MeshBodies.jl")

L = 32
body = MeshBody("/home/marin/Workspace/WaterLily/cube.stl";scale=L/4)
# body = MeshBody("/home/marin/Workspace/CardioVascularFlow.jl/example/3D_flapping/Solid/geom.inp";map,boundary=false,thk=2,scale=1.0)
# struct Tree{T}
#     indices :: Array{}
#     lmax :: Int
#     hyper_rec :: Rect
# end

# outside(x::SVector,bbox::Rect) = !(all(bbox.origin .≤ x) && all(x .≤ bbox.origin+bbox.widths)) # 1.679 ns (0 allocations: 0 bytes)

# split(bbox)

# function search!(tree::Tree{T}, point::SVector{3,T}) where T
#     # are we outside the hyper_rec?
#     outside(point,tree.hyper_rec) && return (dist(point,body.hyper_rec),zero(x),zero(x))

#     # if we are inside, we search children nodes
#     for l in 1:tree.lmax
#         outside(point,tree.children[l].bbox) && (dist(point,body.bbox),zero(x),zero(x))  

#     end
# end