using WaterLily,WaterLilyMeshBodies,StaticArrays,WriteVTK
using FileIO, MeshIO

import WaterLily: @loop
struct ActuatorDisk{T,F<:Function}
    R::T
    map::F
    n::AbstractVector
    function ActuatorDisk(R; map=(x,t)->x)
        new{typeof(R),typeof(map)}(R, map)
    end
end
@fastmath @inline r(x) = √(x[2]^2+x[3]^2)
@fastmath @inline dist(x::SVector, disk::ActuatorDisk) = √(x[1]^2+(r(x)-min(r(x),))^2) - 2.f0
function flux!(Ii::CartesianIndex,disk::ActuatorDisk{T},t) where T
    # map to the location of the disk
    x = disk.map(loc(Ii,T),t); i = Base.last(Ii)
    # find disk normal at this point
    n = 
    return n[i]*WaterLily.kern(clamp(d,-1,1))
end
function force_disk!(flow, t; disk)
    @loop flow.f[Ii] += flux(Ii,disk,t) over Ii in CartesianIndices(flow.u)
end

function make_sphere(;L=32,Re=250,U=1,mem=Array,T=Float32)
    # make the body from the stl mesh
    body = MeshBody("/home/marin/Workspace/CardioVascularFlow.jl/example/aorta.inp";scale=L/2.f0,
                    map=(x,t)->x-SA[L/2.f0,L/2.f0,L/4.f0],boundary=false,half_thk=2.f0,mem)
    # generate sim
    Simulation((L,L,L÷2), (0,0,0), L; ν=U*L/Re, body, mem, T)
end

# make a writer with some attributes to output to the file
vtk_velocity(a::Simulation) = a.flow.u  |> Array;
vtk_pressure(a::Simulation) = a.flow.p  |> Array;
vtk_mu0(a::Simulation)      = a.flow.μ₀ |> Array;
vtk_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, sim_time(a); fastd²=10); a.flow.σ |> Array);
custom_attrib = Dict("u"=>vtk_velocity, "p"=>vtk_pressure, "d"=>vtk_body, "μ₀"=>vtk_mu0)

# save it
mesh_velocity(a) = [SVector(sum(tri,dims=2)/3) for tri in Array(a.velocity)]
wr_mesh = vtkWriter("Aorta_mesh", attrib=Dict("velocity"=>mesh_velocity))

# make the sim
using CUDA
sim = make_sphere(L=64;mem=CuArray)

# make the paraview writer
wr = vtkWriter("Aorta";attrib=custom_attrib)
# duration and write steps
t₀,duration,step = 0.,0.0,0.1

# run the sim
@time for tᵢ in range(t₀,t₀+duration;step)
    # update until time tᵢ in the background
    sim_step!(sim,tᵢ;remeasure=false)
    save!(wr, sim); save!(wr_mesh, sim.body)
    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr); close(wr_mesh)
