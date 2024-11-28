using WaterLily
using MeshIO
using StaticArrays
using WriteVTK
using NearestNeighbors

include("../src/WaterLily/MeshBodies.jl")

function make_sphere(;L=32,Re=250,U=1)
    # make a body
    map(x,t) = x .- SA[L,L,L/2]
    # body = MeshBody("/home/marin/Workspace/WaterLily/cube.stl";map,scale=L/4)
    # body = MeshBody("/home/marin/Workspace/WaterLily/cube.inp";map,scale=L/6)
    # body = MeshBody("/home/marin/Workspace/CalculiX/LIMO_heart/test_vol/geom.inp";map,scale=L/6)
    # @assert all(volume(body) .≈ (L/3)^3) # the volume is exact here!
    body = MeshBody("/home/marin/Workspace/CardioVascularFlow.jl/example/3D_flapping/Solid/geom.inp";map,boundary=false,thk=2,scale=1.0)
    # generate sim
    Simulation((2L,2L,L), (U,0,0), L; ν=U*L/Re, body)
end

bb = AutoBody((x,t)->√sum(abs2,x)-32)

# make a writer with some attributes
velocity(a::Simulation) = a.flow.u |> Array;
pressure(a::Simulation) = a.flow.p |> Array;
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); 
                        a.flow.σ |> Array;)
_vbody(a::Simulation) = a.flow.V |> Array;
mu0(a::Simulation) = a.flow.μ₀ |> Array;
norm(a::Simulation) = (a.flow.f .= 0;
                      @WaterLily.loop a.flow.f[I,:] .= measure(a.body,loc(0,I),0.0)[2] over I ∈ WaterLily.inside(a.flow.p);
                      a.flow.f |> Array;)

custom_attrib = Dict(
    "u" => velocity, "p" => pressure, "n" => norm,
    "d" => _body, "v" => _vbody, "μ₀" => mu0,
)# this

# run a sim and plot the time evolution
sim = make_sphere(L=64)

using BenchmarkTools
body = sim.body;
mesh = body.mesh;
x = SA[20.f0,34.5f0,16.5f0]
t = 0.f0

function test_measure(body,x,t)
    d,n = measure(body.mesh,x,t)
end

# function WaterLily.measure(body::MeshBody,x,t,;fastd²=Inf)
#     # eval d=map(x,t)-x, and n̂
#     # ξ = body.map(x,t)
#     d = zero(eltype(x))
#     n = zero(x)
#     # # #  if we are outside of the bounding box, we can measure approx
#     # bbox = body.bbox;
#     # outside(ξ,bbox) #&& return (dist(ξ,body.bbox),zero(x),zero(x))
#     # d,n = measure(body.mesh,ξ,t)
#     !body.boundary && (d = abs(d)-body.half_thk) # if the mesh is not a boundary, we need to adjust the distance

#     # The velocity depends on the material change of ξ=m(x,t):
#     #   Dm/Dt=0 → ṁ + (dm/dx)ẋ = 0 ∴  ẋ =-(dm/dx)\ṁ
#     J = ForwardDiff.jacobian(x->body.map(x,t), x)
#     dot = ForwardDiff.derivative(t->body.map(x,t), t)
#     return (d,n,-J\dot)
# end

wr = vtkWriter("CalculiX_test"; attrib=custom_attrib)
t₀,duration,step = 0.,10,0.1
write!(wr,sim)
println("tU/L=",round(t₀,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
# @time for tᵢ in range(t₀,t₀+duration;step)
#     # update until time tᵢ in the background
#     sim_step!(sim,tᵢ;remeasure=false)
#     write!(wr,sim)
    
#     # print time step
#     println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
# end
close(wr)