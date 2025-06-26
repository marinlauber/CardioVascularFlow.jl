using WaterLily
using StaticArrays
using WriteVTK
import LinearAlgebra: tr
include("../src/utils.jl")

"""
    Qcriterion2(I::CartesianIndex{3},u)

Q-criterion is a deformation tensor metric to identify vortex cores.
Also see Jeong, J., & Hussain, F., doi:[10.1017/S0022112095000462](https://doi.org/10.1017/S0022112095000462)
"""
function Qcriterion(I::CartesianIndex{3},u)
    J = @SMatrix [WaterLily.∂(i,j,I,u) for i ∈ 1:3, j ∈ 1:3]
    S,Ω = (J+J')/2,(J-J')/2 
    ## -0.5*sum(eigvals(S^2+Ω^2)) # this is also possible, but 2x slower
    0.5*(√(tr(Ω*Ω'))^2-√(tr(S*S'))^2)
end

S(I::CartesianIndex{3},u) = @SMatrix [0.5*(WaterLily.∂(i,j,I,u)+WaterLily.∂(j,i,I,u)) for i ∈ 1:3, j ∈ 1:3]

# seconda invariant viscous stress tensor
function Π₂(I::CartesianIndex{3},u;μ=1)
    # τ = 2*μ*S(I,u) # shear stress tensor from rate of strain
    τ = @SMatrix [μ*WaterLily.∂(i,j,I,u) for i ∈ 1:3, j ∈ 1:3]
    0.5*(tr(τ)^2 - tr(τ^2))
end

function VonMisses(I::CartesianIndex{3},u;μ=1)
    I₂ = Π₂(I,u;μ)
    √(3I₂)
end

function VonMisses2(I,u)
    τ = @SVector [WaterLily.∂(i,j,I,u) for (i,j) in zip((1,2,3),(2,3,1))]
    √(3*sum(abs2,τ))
end

function heart(L=2^5;Re=5e3,mem=Array,U=1,AR=2,T=Float32)
    ## Define simulation size, geometry dimensions, & viscosity
    
    # Motion functions
    ω = T(U/L)
    @fastmath w(t) = 0.1f0*L*cos(ω*t)
    @fastmath h(t) = 0.1f0*L*cos(ω*t)
    @fastmath r(x) = √(x[1]^2+x[2]^2)
    # Build the heart from a mapped ellipsoid and plane
    ellipsoid = AutoBody((x,t)->abs(ellipse(x, SA[L,L/AR]))-1.5f0,
                         (x,t)->x.-SA[-h(t),-h(t),1.5f0*L+w(t)])
    disk = AutoBody((x,t)->√(x[3]^2+(r(x)-min(r(x),L/AR))^2)-1.5f0,
                    (x,t)->x.-SA[-w(t),-w(t),1.85f0*L])
    plane = AutoBody((x,t)->x[3]-1.85f0*L)
    mitral = AutoBody((x,t)->√(x[3]^2+(r(x)-min(r(x),1))^2)-L/6.f0,
                      (x,t)->x.-SA[0,0,1.85f0*L])
    body =  ellipsoid ∩ plane + disk - mitral
    # outer ellipse for visualization
    mask = AutoBody((x,t)->ellipse(x, SA[L,L/AR])-1.5f0,
                     (x,t)->x.-SA[-h(t),-h(t),1.5f0*L+w(t)]) 

    # Return initialized simulation
    Simulation((Int(0.75*L),Int(0.75*L),3*L),(0f0,0f0,0f0),L;U,ν=U*L/Re,body,mem,T=T),mask
end

using CUDA
CUDA.allowscalar(false)
# Define geometry and motion on GPU    
sim,mask = heart(3*2^5;mem=CUDA.CuArray);

# make a writer with some attributes
vtk_velocity(a::Simulation) = a.flow.u |> Array;
vtk_pressure(a::Simulation) = a.flow.p |> Array;
vtk_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); 
                           a.flow.σ |> Array;)
vtk_vorticity(a::Simulation) = (@inside a.flow.σ[I] = WaterLily.curl(3,I,a.flow.u)*a.L/a.U;
                                a.flow.σ |> Array;)
vtk_vbody(a::Simulation) = a.flow.V |> Array;
vtk_mu0(a::Simulation) = a.flow.μ₀ |> Array;
vtk_lamda(a::Simulation) = (measure_sdf!(a.flow.σ, mask, WaterLily.time(a.flow));
                            @inside a.flow.σ[I] = ifelse(a.flow.σ[I]<0,WaterLily.λ₂(I,a.flow.u),0);
                            a.flow.σ |> Array;)
vtk_fake_outside(a::Simulation) = (measure_sdf!(a.flow.σ, mask, WaterLily.time(a.flow));
                                   a.flow.σ |> Array;)
vtk_vbody(a::Simulation) = a.flow.V |> Array;
vtk_Q(a) = (@inside a.flow.σ[I] = Qcriterion(I,a.flow.u); a.flow.σ |> Array)


custom_attrib = Dict( "u" => vtk_velocity, "p" => vtk_pressure,  "λ₂" => vtk_lamda, "Q" => vtk_Q)# this


wr = vtkWriter("test"; attrib=custom_attrib)

# Loop in time
foreach(1:200) do frame
    @show frame
    sim_step!(sim,sim_time(sim)+0.05,verbose=true);
    write!(wr, sim)
end
close(wr)
