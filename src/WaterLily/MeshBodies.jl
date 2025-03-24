using FileIO, MeshIO
using GeometryBasics
using WaterLily: AbstractBody,@loop,measure
using StaticArrays
using ForwardDiff
# using NearestNeighbors

struct MeshBody{T,F<:Function} <: AbstractBody
    mesh  :: GeometryBasics.Mesh
    map   :: F
    bbox  :: Rect
    scale :: T
    half_thk::T #half thickness
    boundary::Bool
end
function MeshBody(fname;map=(x,t)->x,scale=1.0,boundary=true,thk=0f0,T=Float32)
    tmp = endswith(fname,".inp") ? load_inp(fname) : load(fname)
    points = GeometryBasics.Point.(tmp.position*scale) # can we specify types?
    mesh = GeometryBasics.Mesh(points,GeometryBasics.faces(tmp))
    bbox = Rect(mesh.position)
    bbox = Rect(bbox.origin.-max(4,thk),bbox.widths.+max(8,2thk))
    MeshBody(mesh,map,bbox,T(scale),T(thk/2),boundary)
end
Base.copy(b::MeshBody) = MeshBody(GeometryBasics.Mesh(b.mesh.position,GeometryBasics.faces(b.mesh)),
                                  b.map,Rect(b.bbox),b,scale,b.half_thk,b.boundary)

function load_inp(fname; facetype=GLTriangleFace, pointtype=Point3f)
    #INP file format
    @assert endswith(fname,".inp") "file type not supported"
    fs = open(fname)

    points = pointtype[]
    faces = facetype[]
    node_idx = Int[]

    # read the first 3 lines if there is the "*heading" keyword
    line = readline(fs)
    contains(line,"*heading") && (line = readline(fs))
    BlockType = contains(line,"*NODE") ? Val{:NodeBlock}() : Val{:DataBlock}()
    
    # read the file
    while !eof(fs)
        line = readline(fs)
        BlockType, line = parse_blocktype!(BlockType, fs, line)
        if BlockType == Val{:NodeBlock}()
            push!(node_idx, parse(Int,split(line,",")[1])) # keep track of the node index of the inp file
            push!(points, pointtype(parse.(eltype(pointtype),split(line,",")[2:4])))
        elseif BlockType == Val{:ElementBlock}()
            nodes = parse.(Int,split(line,",")[2:end])
            push!(faces, TriangleFace{Int}(facetype([findfirst(==(node),node_idx) for node in nodes])...)) # parse the face
        else
            continue
        end
    end
    return Mesh(points, faces); close(fs);
end
function parse_blocktype!(block, io, line)
    contains(line,"*NODE") && return block=Val{:NodeBlock}(),readline(io)
    contains(line,"*ELEMENT") && return block=Val{:ElementBlock}(),readline(io)
    contains(line,"*ELSET, ELSET=") && return block=Val{:SkipBlock}(),readline(io)
    return block, line
end

"""
    locate(p,tri)

Find the closest point `x` on the triangle `tri` to the point `p`.
"""
function locate(tri::GeometryBasics.Ngon{3},p::SVector{T}) where T #5.327 ns (0 allocations: 0 bytes)
    # is point `a` closest?
    a, b, c = tri.points
    ab = b.-a
    ac = c.-a
    ap = p.-a
    d1 = sum(ab.*ap)
    d2 = sum(ac.*ap)
    # is point `a` closest?
    if ((d1 ≤ 0) && (d2 ≤ 0))
        return a
    end
    # is point `b` closest?
    bp = p.-b
    d3 = sum(ab.*bp)
    d4 = sum(ac.*bp)
    if ((d3 ≥ 0) && (d4 ≤ d3))
        return b
    end
    # is point `c` closest?
    cp = p.-c
    d5 = sum(ab.*cp)
    d6 = sum(ac.*cp)
    if ((d6 ≥ 0) && (d5 ≤ d6))
        return c
    end
    # is segment 'ab' closest?
    vc = d1*d4 - d3*d2
    if ((vc ≤ 0) && (d1 ≥ 0) && (d3 ≤ 0))
        x =  a .+ ab.*d1 ./ (d1 - d3)
        return x
    end
    #  is segment 'ac' closest?
    vb = d5*d2 - d1*d6
    if ((vb ≤ 0) && (d2 ≥ 0) && (d6 ≤ 0))
        x =  a .+ ac.*d2 ./ (d2 - d6)
        return x
    end
    # is segment 'bc' closest?
    va = d3*d6 - d5*d4
    if ((va ≤ 0) && (d4 ≥ d3) && (d5 ≥ d6))
        x =  b .+ (c .- b) .* (d4 - d3) ./ ((d4 - d3) + (d5 - d6))
        return x
    end
    # closest is interior to `abc`
    denom = one(T) / (va + vb + vc)
    v= vb*denom
    w = vc*denom
    x = a .+ ab .* v .+ ac .* w
    return x
end

# inside bbox or not
outside(x::SVector,bbox::Rect) = !(all(bbox.origin .≤ x) && all(x .≤ bbox.origin+bbox.widths)) # 1.679 ns (0 allocations: 0 bytes)
rect = Rect(0,0,0,1,1,1) # origin and widths
# @assert !inside(SA[0.5,1,2.5],rect) && inside(SA[0.5,0.5,0.5],rect)
# 
# distance to box center
dist(x::SVector,bbox::Rect) = √sum(abs2,x.-bbox.origin-0.5bbox.widths) # 1.707 ns (0 allocations: 0 bytes)
rect = Rect(0,0,0,1,1,1) # origin and widths
# @assert dist(SA[1.0,1.0,1.0],rect) == √0.75
# @assert dist(SA[1.5,1.0,1.0],rect) == √1.5
# @assert dist(SA[1.5,1.5,1.0],rect) == √2.25
# @assert dist(SA[1.5,1.5,1.5],rect) == √3.0

WaterLily.sdf(body::MeshBody,x,t;kwargs...) = (d=measure(body.mesh,body.map(x,t),t;kwargs...)[1];
                                               !body.boundary && (d = abs(d)-body.half_thk); return d) # if the mesh is not a boundary, we need to adjust the distance

function WaterLily.measure(body::MeshBody,x,t,;kwargs...)
    # eval d=map(x,t)-x, and n̂
    ξ = body.map(x,t)
    #  if we are outside of the bouding box, we can measure approx
    outside(ξ,body.bbox) && return (dist(ξ,body.bbox),zero(x),zero(x))
    d,n = measure(body.mesh,ξ,t)
    !body.boundary && (d = abs(d)-body.half_thk) # if the mesh is not a boundary, we need to adjust the distance

    # The velocity depends on the material change of ξ=m(x,t):
    #   Dm/Dt=0 → ṁ + (dm/dx)ẋ = 0 ∴  ẋ =-(dm/dx)\ṁ
    J = ForwardDiff.jacobian(x->body.map(x,t), x)
    dot = ForwardDiff.derivative(t->body.map(x,t), t)
    return (d,n,-J\dot)
end
# function measure_nn(body::MeshBody,x,t;kwargs...)
#     # eval d=map(x,t)-x, and n̂
#     ξ = body.map(x,t)
#     # measure with the KDTree
#     u,_ = nn(@views(body.tree), ξ)
#     v = x-locate(body.mesh[u],x)
#     n = normal(body.mesh[u])
#     d = sign(sum(v.*n))*√sum(abs2,v) # signed Euclidian distance
#     !body.boundary && (d = abs(d)-body.half_thk) # if the mesh is not a boundary, we need to adjust the distance

#     # The velocity depends on the material change of ξ=m(x,t):
#     #   Dm/Dt=0 → ṁ + (dm/dx)ẋ = 0 ∴  ẋ =-(dm/dx)\ṁ
#     J = ForwardDiff.jacobian(x->body.map(x,t), x)
#     dot = ForwardDiff.derivative(t->body.map(x,t), t)
#     return (d,n,-J\dot)
# end

using LinearAlgebra: cross
"""
    normal(tri::GeometryBasics.Ngon{3})

Return the normal vector to the triangle `tri`.
"""
function normal(tri::GeometryBasics.Ngon{3}) #3.039 ns (0 allocations: 0 bytes)
    n = cross(SVector(tri.points[2]-tri.points[1]),SVector(tri.points[3]-tri.points[1]))
    n/√sum(abs2,n)
end
@inbounds @inline d²_fast(tri::GeometryBasics.Ngon{3},x) = sum(abs2,x-center(tri)) #1.825 ns (0 allocations: 0 bytes)
@inbounds @inline d²(tri::GeometryBasics.Ngon{3},x) = sum(abs2,x-locate(tri,x)) #4.425 ns (0 allocations: 0 bytes)
@inbounds @inline center(tri::GeometryBasics.Ngon{3}) = SVector(sum(tri.points;dims=1)/3.f0...) #1.696 ns (0 allocations: 0 bytes)
@inbounds @inline area(tri::GeometryBasics.Ngon{3}) = 0.5*√sum(abs2,cross(SVector(tri.points[2]-tri.points[1]),SVector(tri.points[3]-tri.points[1]))) #1.784 ns (0 allocations: 0 bytes)
"""
    measure(mesh::GeometryBasics.Mesh,x,t;kwargs...)

Measure the distance `d` and normal `n` to the mesh at point `x` and time `t`.
"""
function WaterLily.measure(mesh::GeometryBasics.Mesh,x::SVector{T},t;kwargs...) where T
    u=1; a=b=d²_fast(mesh[1],x) # fast method
    for I in 2:length(mesh)
        b = d²_fast(mesh[I],x)
        b<a && (a=b; u=I) # Replace current best
    end
    n = normal(mesh[u])
    v = x-locate(mesh[u],x) # signed Euclidian distance
    d = sign(sum(v.*n))*√sum(abs2,v) 
    return d,n
end # 120.029 ns (0 allocations: 0 bytes)d # 4.266 μs (0 allocations: 0 bytes)

# use divergence theorem to calculate volume of surface mesh
# F⃗⋅k⃗ = -⨕pn⃗⋅k⃗ dS = ∮(C-ρgz)n⃗⋅k⃗ dS = ∫∇⋅(C-ρgzk⃗)dV = ρg∫∂/∂z(ρgzk⃗)dV = ρg∫dV = ρgV #
function volume(body::MeshBody)
    vol = zeros(3)
    for tri in body.mesh
        # this is the hydrostatic pressure contribution of every triangle in the 3 axis
        vol .+= center(tri) .* normal(tri) .* area(tri)
    end
    return vol
end

import WaterLily: interp
function forces(body::MeshBody,flow::Flow)
    f = zeros((3,len(body.mesh)))
    for (i,tri) in enumerate(body.mesh)
        f[:,i] .= normal(tri) .* area(tri) * WaterLily.interp(center(tri) .* normal(tri), flow.p)
    end
    return f
end