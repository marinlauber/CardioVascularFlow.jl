using WaterLily, StaticArrays, WaterLilyMeshBodies, GeometryBasics
using GeometryBasics: Mesh, Point3f, TriangleFace, GLTriangleFace
using LinearAlgebra: dot

# limit the time step
WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.25)

function make_sim(;L=64,Re=5000,U=1,mem=Array,T=Float32)
    # extract the good stuff
    data = open(joinpath(@__DIR__, "endocardium_motion.txt"), "r") do file
        mapreduce(L->parse.(T, split(L)), hcat, readlines(file)[10:end])
    end
    # remove the reference height
    zmin,zmax = extrema(data[3,:]); scale = T(zmax-zmin)
    # maximum radius
    R_max = max(diff([extrema(data[1,:])...])[1] , diff([extrema(data[2,:])...])[1])/2.f0/scale*L

    # remove 3 first columns since these are the ref config, make it unit-length and scale
    motion = [data[1+i*3:3+i*3,:]./scale*L for i in 1:size(data)[1]÷3-1]
    # the top surface moves up and down, we must fix it
    motion = [m.+SA{T}[1L/2.f0,1L/2.f0,5L÷4-maximum(m[3,:])] for m in motion]

    # read it from the connecivity file
    connectivity = open(joinpath(@__DIR__, "connectivity.txt"), "r") do f
        mapreduce(L->parse.(Int, split(L)), hcat, readlines(f))
    end

    # make points list from first snapshot, we use the last config such that
    # measuring first sets the correct velocity BC
    points = [Point3f(pnt) for pnt in eachcol(motion[end])]

    # initial mesh, create triangulation
    # connectivity must be incremented since we have 1 indexing here
    faces = [TriangleFace{Int}(GLTriangleFace((reverse(tri).+1...))) for tri in eachcol(connectivity)]

    # generate a mesh from this
    mesh = GeometryBasics.Mesh(points, faces)
    shell = MeshBody(mesh;scale=1.f0,half_thk=2f0,boundary=false,mem=mem)

    # add a top surface
    R₁,R₂ = T(1.25*R_max),T(L/9)
    println("R₁=",R₁,"($(typeof(R₁))), R₂=",R₂,"($(typeof(R₂)))")
    map(x,t) = x.-SA[L/2.f0,L/2.f0,5.f0L÷4.f0]
    cap =  AutoBody((x,t)->√(x[3]^2+(√(x[1]^2+x[2]^2)-min(√(x[1]^2+x[2]^2),R₁))^2)-2.0f0, map)
    cap -= AutoBody((x,t)->√(x[3]^2+(√(x[1]^2+x[2]^2)-min(√(x[1]^2+x[2]^2),R₂))^2)-4.0f0, map)

    # make the body
    body = shell + cap

    # make motion usable for interpolate, and split into cylindrical components
    # about the vertical axis through (L/2,L/2): the displacement of each vertex
    # relative to the first snapshot is projected onto the local cylindrical basis
    # (êr,êθ,êz), giving a vertical (z), an in-plane radial (r) and a twisting (θ)
    # contribution that sum back to Δp = motion[i] - motion[1].
    N_snaps = length(motion)
    N_elems = length(faces)
    ref_motion = Array{SMatrix{3,3,T,9}}(undef, N_snaps, N_elems)
    z_motion = similar(ref_motion)
    r_motion = similar(ref_motion)
    θ_motion = similar(ref_motion)
    cx, cy, êz = T(L/2.f0), T(L/2.f0), SA{T}[0,0,1]
    # decompose a vertex displacement Δ at reference position p₀ → (z,r,θ) vectors
    cyl(p₀, Δ) = begin
        dx, dy = p₀[1]-cx, p₀[2]-cy
        ρ = √(dx^2+dy^2)
        êr = ρ > T(1e-6) ? SA{T}[dx/ρ, dy/ρ, 0] : zero(SVector{3,T})
        êθ = ρ > T(1e-6) ? SA{T}[-dy/ρ, dx/ρ, 0] : zero(SVector{3,T})
        (dot(Δ,êz)*êz, dot(Δ,êr)*êr, dot(Δ,êθ)*êθ)
    end
    for i ∈ 1:N_snaps
        pts  = [Point3f(pnt) for pnt in eachcol(motion[i])]
        msh  = GeometryBasics.Mesh(pts, faces)
        for j in 1:N_elems
            cur = SMatrix{3,3,T}(hcat(msh[j]...))
            ref_motion[i, j] = cur
            ref = ref_motion[1, j]  # first snapshot is the reference config
            c = ntuple(k -> cyl(ref[:,k], cur[:,k]-ref[:,k]), 3)  # per-vertex (z,r,θ)
            z_motion[i, j] = hcat(c[1][1], c[2][1], c[3][1])
            r_motion[i, j] = hcat(c[1][2], c[2][2], c[3][2])
            θ_motion[i, j] = hcat(c[1][3], c[2][3], c[3][3])
        end
    end
    ref_motion = mem(ref_motion)

    # make a simulation
    Simulation((L,L,2L),(0,0,0),L;U,ν=U*L/Re,body,mem,T), ref_motion, shell,
        (z=mem(z_motion), r=mem(r_motion), θ=mem(θ_motion))
end

using CUDA
# make simulation
sim,ref_motion,mesh_body,cyl_motion = make_sim(L=96;mem=CuArray)
@assert all(isapprox.(ref_motion, ref_motion[1:1, :] .+ (cyl_motion.z .+ cyl_motion.r .+ cyl_motion.θ); atol=1f-4))

# build a motion by combining z,r,θ motions
R=false
Z=true
θ=true
motion_data = ref_motion[1:1, :] .+ (Int(Z).*cyl_motion.z .+ Int(R).*cyl_motion.r .+ Int(θ).*cyl_motion.θ)

# make a motion interpolation from the snapshot
times = collect(range(0.f0, sim.L, length=size(motion_data, 1)))
motion = MotionInterpolation(motion_data, times; periodic=true)

# flow writer
vtk_velocity(a::AbstractSimulation) = a.flow.u  |> Array;
vtk_μ₀(a::AbstractSimulation)       = a.flow.μ₀ |> Array;
vtk_pressure(a::AbstractSimulation) = a.flow.p  |> Array;
vtk_vbody(a::AbstractSimulation)    = a.flow.V  |> Array;
vtk_ω(a::AbstractSimulation)        = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u);     a.flow.σ |> Array);
vtk_lambda(a::AbstractSimulation)   = (@inside a.flow.σ[I] = WaterLily.λ₂(I,a.flow.u);         a.flow.σ |> Array);
vtk_body(a::AbstractSimulation)     = (measure_sdf!(a.flow.σ, a.body, sim_time(a); fastd²=10); a.flow.σ |> Array);

# writer for the flow
new_attrib = Dict("u"=>vtk_velocity, "μ₀"=>vtk_μ₀,
                  "V"=>vtk_vbody, "p"=>vtk_pressure,
                  "d"=>vtk_body, "λ₂"=>vtk_lambda, "ω₃"=>vtk_ω)
wr = vtkWriter("endocardium_flow"; attrib=new_attrib)

# writer for the structure
mesh_velocity(a) = [SVector(sum(tri,dims=2)/3) for tri in Array(a.velocity)]
wr_msh = vtkWriter("endocardium_mesh", attrib=Dict("velocity"=>mesh_velocity))

# Run for a few cycles
for tᵢ in range(0, 1, step=0.02)
    while sim_time(sim) < tᵢ
        push!(sim.flow.Δt, 0.5)
    end
    @show tᵢ
    # get the position an mesh velocity at this time
    sim.body = interpolate!(sim.body, motion, WaterLily.time(sim))
    # measure the new body
    measure!(sim)
    # save both the mesh and the flow
    save!(wr_msh, sim.body, sim_time(sim)); save!(wr, sim)
end
close(wr_msh); close(wr)