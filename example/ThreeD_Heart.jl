using WaterLily
using StaticArrays
using WriteVTK
include("../../Tutorials-WaterLily/src/ThreeD_plots.jl")
include("../src/utils.jl")

function heart(L=2^5;Re=5e2,mem=Array,U=1,AR=2)
    # Define simulation size, geometry dimensions, & viscosity
    
    # Motion functions
    ω = 2U/L
    @fastmath w(t) = 0.1*L*cos(ω*t)
    @fastmath h(t) = 0.1*L*cos(ω*t)
    @fastmath r(x) = √(x[1]^2+x[2]^2)
    # Build the heart from a mapped ellipsoid and plane
    ellipsoid = AutoBody((x,t)->abs(ellipse(x, SA[L,L/AR]))-1.5,
                         (x,t)->x.-SA[-h(t),-h(t),1.5L+w(t)])
    disk = AutoBody((x,t)->√(x[3]^2+(r(x)-min(r(x),L/AR))^2)-1.5,
                    (x,t)->x.-SA[-w(t),-w(t),1.85L])
    plane = AutoBody((x,t)->x[3]-1.85L)
    mitral = AutoBody((x,t)->√(x[3]^2+(r(x)-min(r(x),1))^2)-5,
                      (x,t)->x.-SA[0,0,1.85L])
    body =  ellipsoid ∩ plane + disk - mitral

    # Return initialized simulation
    Simulation((Int(0.75*L),Int(0.75*L),3L),(0,0,U),L;U,ν=U*L/Re,body,mem,T=Float32)
end

function mirrorto!(a,b)
    n = size(b,1)
    a[reverse(1:n),reverse(1:n),:].=b
    a[reverse(n+1:2n),1:n,:].=a[1:n,1:n,:]
    a[:,reverse(n+1:2n),:].=a[:,1:n,:]
    return a
end
Makie.inline!(false)
# using CUDA
# CUDA.allowscalar(false)
begin
    # Define geometry and motion on GPU
    sim = heart()#mem=CUDA.CuArray);
    # sim_step!(sim,sim_time(sim)+0.05);

    # Create CPU buffer arrays for geometry flow viz 
    a = sim.flow.σ
    d = similar(a,size(inside(a))) |> Array; # one quadrant
    md = similar(d,(2,2,1).*size(d))  # hold mirrored data

    # Set up geometry viz
    geom = geom!(md,d,sim) |> Observable;
    fig, _, _ = GLMakie.mesh(geom, alpha=0.1, color=:red)

    #Set up flow viz
    ω = ω!(md,d,sim) |> Observable;
    volume!(ω, algorithm=:mip, colormap=:algae, colorrange=(1,10))
    fig
end

# WaterLily.CFL(a::Flow) = WaterLily.CFL(a;Δt_max=0.2)

wr = vtkWriter("heart"; attrib=custom_attrib)

# Loop in time
# record(fig,"jelly.mp4",1:200) do frame
foreach(1:100) do frame
    @show frame
    sim_step!(sim,sim_time(sim)+0.05,verbose=true);
    write!(wr, sim)
    geom[] = geom!(md,d,sim);
    ω[] = ω!(md,d,sim);
end
close(wr)