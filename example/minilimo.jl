using WaterLily, StaticArrays, WaterLilyMeshBodies, JLD2
using GeometryBasics, CUDA, WriteVTK

function make_sim(;L=64,Re=5000,U=1,mem=CuArray,T=Float32)
    # extract the good stuff
    data = jldopen(joinpath("/home/marin/Workspace/FerriteShells.jl/", "limo_motion.jld2"), "r")
    coordinates = data["coordinates"]
    connectivity = data["connectivity"]

    # remove the reference height
    xmin,xmax = extrema(coordinates[1][1,:]); scale = xmax-xmin
    coordinates ./= scale # make it unit length

    # make points list from first snapshot, we use the last config such that
    # measuring first sets the correct velocity BC
    points = [Point3f(pnt) for pnt in eachcol(coordinates[1])]

    # initial mesh, create triangulation from the Quad arrays
    faces = [QuadFace(tri) for tri in eachcol(connectivity)]

    # Build motion as a 2D Array: [N_snapshots × N_elements]
    faces_tri = GeometryBasics.decompose(GLTriangleFace, faces)
    N_snaps = length(coordinates)
    N_elems = length(faces_tri)
    motion_data = Array{SMatrix{3,3,T,9}}(undef, N_snaps, N_elems)
    for i ∈ 1:N_snaps
        pts  = [Point3f(pnt * L) for pnt in eachcol(coordinates[i])]
        msh  = GeometryBasics.Mesh(pts, faces_tri)
        for j in 1:N_elems
            motion_data[i, j] = SMatrix{3,3,T}(hcat(msh[j]...))
        end
    end
    motion_data = mem(motion_data)

    # map to center
    function map(x,t)
        y = x .- SA[L, 5L/8.f0, L/2.f0]
        return SA[y[1], y[2], abs(y[3])]
    end

    # generate a mesh from this
    mesh = GeometryBasics.Mesh(points, faces)
    shell = MeshBody(mesh;map=map,scale=T(L),half_thk=2f0,boundary=false,mem=mem)

    # add a top surface
    idx = findall(≈(xmin/scale), coordinates[1][1,:]) # top edge node indices
    ymin,ymax = extrema(coordinates[1][2,idx])
    zmin,zmax = extrema(coordinates[1][3,idx])
    R₁,R₂ = T((ymax-ymin)*L/2), T((zmax-zmin)*L/2)
    map_cap(x,t) = x.-SA[L/2.f0.-2.f0,5L/8.f0,L/2.f0]
    cap = AutoBody((x,t)->√(x[1]^2+(√(x[2]^2+(x[3]*R₁/R₂/2)^2)-min(√(x[2]^2+(x[3]*R₁/R₂/2)^2),R₁))^2)-2.0f0, map_cap)
    # the holes
    Rᵥ = T(0.9R₂)
    cap -= AutoBody((x,t)->√(x[1]^2+(√(x[2]^2+x[3]^2)-min(√(x[2]^2+x[3]^2),Rᵥ))^2)-4.0f0, map_cap)

    # make the body
    body = shell + cap

    # make a simulation
    sim = Simulation((2L,Int(5L÷4),L),(0,0,0),L;U,ν=U*L/Re,body,mem,T)

    return sim, motion_data, coordinates[1]
end

# params
T = 1.0
sim, motion_data, c = make_sim(;L=128)
times = collect(range(0.f0, sim.L, length=size(motion_data, 1)))

motion = MotionInterpolation(motion_data, times; periodic=true)

# flow writer
vtk_velocity(a::AbstractSimulation) = a.flow.u |> Array;
vtk_μ₀(a::AbstractSimulation) = a.flow.μ₀ |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p |> Array;
vtk_ω(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u); a.flow.σ |> Array);
vtk_vbody(a::AbstractSimulation) = a.flow.V |> Array;
vtk_lambda(a::AbstractSimulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I,a.flow.u); a.flow.σ |> Array);
vtk_body(a::AbstractSimulation) = (measure_sdf!(a.flow.σ, a.body, sim_time(a); fastd²=10); a.flow.σ |> Array);
new_attrib = Dict("u"=>vtk_velocity, "μ₀"=>vtk_μ₀, "V"=>vtk_vbody, "p"=>vtk_pressure, "d"=>vtk_body, "λ₂"=>vtk_lambda, "ω₃"=>vtk_ω)
wr = vtkWriter("limo_flow"; attrib=new_attrib)

# save it
mesh_velocity(a) = [SVector(sum(tri,dims=2)/3) for tri in Array(a.velocity)]
wr_msh = vtkWriter("limo_mesh", attrib=Dict("velocity"=>mesh_velocity))

# initial state
save!(wr_msh, sim.body, 0.0)
save!(wr, sim)

# interpolate and measure
for tᵢ in range(0, 1, step=0.02)
    while sim_time(sim) < tᵢ
        push!(sim.flow.Δt, 0.5)
    end
    @show tᵢ
    sim.body = interpolate!(sim.body, motion, WaterLily.time(sim))
    measure!(sim)
    save!(wr_msh, sim.body, sim_time(sim))
    save!(wr, sim)
end
close(wr_msh); close(wr)
