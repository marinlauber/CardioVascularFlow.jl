using StaticArrays
using WaterLily
using BenchmarkTools
using StaticArrays
using NearestNeighbors
using Plots

# include("MeshBodies.jl")

# L = 32
# body = MeshBody("/home/marin/Workspace/CardioVascularFlow.jl/example/sphere.inp";scale=L)

# 
# bbox have C (center) and R (width)
struct Bbox{S<:SVector}
    C::S
    R::S
    leaf::Bool
    indices
    function Bbox(C::S,R::S,leaf::Bool=false,indices=[0]) where S
        new{S}(C,R,leaf,indices)
    end
end
box = Bbox(SA[0.0,0.0,0.0],SA[0.5,0.5-eps(),0.5-eps()])
# split the box in the longest side
split_w(width::SVector{N},j) where N = SA[ntuple(i -> i==j ? width[i]/2 : width[i], N)...]
function Base.split(b::Bbox)
    # split the longest side
    w = split_w(b.R, argmax(b.R))
    return Bbox(b.C-(b.R-w),w), Bbox(b.C+(b.R-w),w)
end
left,right = split(box)
@btime split_w($box.R,$(argmax(box.R))) # 2.266 ns , 0 bytes
# @assert all(left.C.≈SA[-0.25,0.,0.]) && all(left.R.≈SA[0.25,0.5,0.5])
# @assert all(right.C.≈SA[0.25,0.,0.]) && all(right.R.≈SA[0.25,0.5,0.5])

# check if a point is outside the box
WaterLily.inside(x,b::Union{Tree,Bbox}) = (all(b.C-1.2b.R .≤ x) && all(x .≤ b.C+1.2b.R))
point = [0.,0.,0.]
b = Bbox(SA[0.0,0.0,0.0],SA[0.5,0.5,0.5])
@assert inside(point,b) # 4.985 ns , 0 bytes
@btime inside($point,$b)
# @assert !inside([0.5+eps(),0.5,0.5],b) && !inside([-0.5-eps(),0.5,0.5],b)

function bounding_box(points::AbstractArray{T,2};δ=1e-2) where T
    vmax = SVector(ntuple(i->typemin(T),size(points,1)))
    vmin = SVector(ntuple(i->typemax(T),size(points,1)))
    for p in eachcol(points)
        vmin = min.(p, vmin)
        vmax = max.(p, vmax)
    end
    o = (vmin + vmax)/2
    r = (vmax - vmin)/2 .+ δ # make it a bit bigger
    return Bbox(SA[o...],SA[r...],false)
end
# some points
points = rand(3,10)
box = bounding_box(points)
left, right = split(box)

# check if points are inside the box
# Base.filter(pts,b::Bbox) = filter(x->inside(x,b),eachcol(pts))
Base.filter(pts,b::Bbox) = findall(x->inside(x,b),eachcol(pts))
p_left = filter(points,left)
p_right = filter(points,right)
# @assert length(p_left) + length(p_right) == size(points,2)
# @assert !(p_left in p_right) # check that there is no point in both boxes

# function tree(points::AbstractArray{T,2}; Nmax=100) where T
#     b = bounding_box(points,δ=0.1)
#     boxes,idx = Bbox[b],Int[1]
#     # split the box in the longest direction
#     left, right = split(b)
#     # find which box each point belongs to
#     points_left = filter(points, left)
#     points_right = filter(points, right) # we can avoid this actually
#     # downsample the boxes
#     tree!(boxes, idx, 2, @views(points[:,points_left]), left, Nmax)
#     tree!(boxes, idx, 3, @views(points[:,points_right]), right, Nmax)
#     #  sort the boxes by index and return
#     return tuple(map(i->boxes[i],sortperm(idx))...)
# end

function tree!(list, idx, Is, points::AbstractArray{T,2}, b::Bbox, Nmax) where T
    push!(idx,Is) # add the index location
    # if we are on a leaf, we push a leaf box on the list
    if (length(points)<Nmax)
        push!(list, Bbox(b.C,1.2*b.R,true,parentindices(points)[2]))
        return nothing
    end
    push!(list, b);
    # we are not on a leaf node, we can downsample
    left, right = split(b)
    # find which box each point belongs to
    points_left = filter(points, left)
    points_right = filter(points, right)
    # downsample the boxes
    tree!(list, idx, Is*2, @views(points[:,points_left]), left, Nmax)
    tree!(list, idx, Is*2+1, @views(points[:,points_right]), right, Nmax)
    return nothing
end

function measure(tree::Tree, x::SVector)
    # are we in the bbox of the geom?
    !inside(x,tree) && return √sum(abs2,max.(0,abs.(x-tree.C)-tree.R))
    # if we are inside the bbox, which sub-box are we in?
    inleft = inside(x,tree[2]) # if we are not in the left box, we are in the right one
    return _measure(tree, inleft ? 2 : 3, x)
end

function _measure(tree::Tree,Is,x::SVector)
    b = tree[Is]
    # if we are on a leaf, we measure the points inside
    b.leaf && return b.indices[1]
    # if not, we check which box we are in
    inleft = inside(x,tree[2Is]) # if we are not in the left box, we are in the right one
    return _measure(tree, inleft ? 2Is : 2Is+1, x)
end

mutable struct Tree
    boxes::NTuple
    points::AbstractArray
    C::SVector # these point to boxes[1]
    R::SVector # these point to boxes[1]
    function Tree(points::AbstractArray{T,2}; Nmax=50) where T
        b = bounding_box(points,δ=0.1)
        boxes,idx = Bbox[b],Int[1]
        # split the box in the longest direction
        left, right = split(b)
        # find which box each point belongs to
        points_left = filter(points, left)
        points_right = filter(points, right) # we can avoid this actually
        # downsample the boxes
        tree!(boxes, idx, 2, @views(points[:,points_left]), left, Nmax)
        tree!(boxes, idx, 3, @views(points[:,points_right]), right, Nmax)
        #  sort the boxes by index and return
        boxes = tuple(map(i->boxes[i],sortperm(idx))...)
        new(boxes,points,boxes[1].C,boxes[1].R)
    end
end
Base.getindex(t::Tree, i::Int64) = t.boxes[i]


# make a simple tree
points = rand(2,100)
points[:,1] .= [-1,-1]
# try a circle
points[1,:] .= cos.(range(0,stop=2π,length=100))
points[2,:] .= sin.(range(0,stop=2π,length=100))

tree = Tree(points; Nmax=50)

d = rand(1:size(tree.boxes,1))
let
    plot(title="KDTree - see $d",dpi=300)
    for (i,b) in enumerate(tree.boxes)
        p1 = b.C - b.R
        p2 = b.C + SA[b.R[1],-b.R[2]]
        p3 = b.C + b.R
        p4 = b.C - SA[b.R[1],-b.R[2]]
        plot!([p1[1],p2[1],p3[1],p4[1],p1[1]],[p1[2],p2[2],p3[2],p4[2],p1[2]],
               color=:black, lw=2,fill=false,alpha=0.3,label=:none)
        c = ifelse(b.leaf,:forestgreen,:brown3)
        annotate!(b.C[1], b.C[2], ("$i", 16, c, :center))
    end
    plot!(points[1,:],points[2,:],seriestype=:scatter,label=:none)
    # select a box and see what's inside
    idx = tree[d].indices
    plot!(points[1,idx],points[2,idx],seriestype=:scatter,label=:none)
end