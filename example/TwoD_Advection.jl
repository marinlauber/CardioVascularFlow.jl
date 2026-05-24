using WaterLily,StaticArrays,Plots
import WaterLily: @loop
using InterfaceAdvection
import InterfaceAdvection: MPCFL

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

# second invariant viscous stress tensor
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

# Two-phase simulation type
mutable struct TwoPhaseSimulation <: AbstractSimulation
    sim :: Simulation
    intf :: cVOF
    function TwoPhaseSimulation(dims, args...; ν=0., perdir=(), λμ=1e-2, λρ=1e-3, η=0,
                                InterfaceSDF=(x)->-5-x[1], mem=Array, kwargs...)
        # WaterLily simulation
        sim = Simulation(dims, args...; ν, perdir, mem, kwargs...)
        # cVOF interface
        intf = cVOF(dims;arr=mem,T=eltype(sim.flow.u),InterfaceSDF,μ=ν,λμ,λρ,η,perdir)
        # print info
        println("μ: $(intf.μ), λρ: $(intf.λρ)")
        # adjust time step
        sim.flow.Δt[end] = min(sim.flow.Δt[end],MPCFL(sim.flow,intf))
        new(sim,intf)
    end
end
# overload properties
Base.getproperty(f::TwoPhaseSimulation, s::Symbol) = s in propertynames(f) ? getfield(f, s) : getfield(f.sim, s)
Base.setproperty!(f::TwoPhaseSimulation, s::Symbol, x) = s in propertynames(f) ? setproperty!(f,s,x) : setproperty!(f.sim,s,x)

import InterfaceAdvection: advectVOF!
function passive_advect!(a::TwoPhaseSimulation;production=nothing)
    u⁰,u,Δt,T = a.flow.u⁰,a.flow.u,a.flow.Δt[end-1],eltype(a.flow.u)
    a.intf.fᶠ .= a.intf.f # set initial condition
    production!(production,u⁰,u,a.intf)
    advectVOF!(a.intf.f,a.intf.fᶠ,a.intf.α,a.intf.n̂,u⁰,u,Δt,a.intf.c̄,a.intf.ρuf,zero(T);perdir=a.flow.perdir)
    # add production term
    a.intf.fᶠ .= (a.intf.fᶠ+a.intf.f)/2 # mid-point for predictor-corrector
    production!(production,u,u,a.intf)
    advectVOF!(a.intf.f,a.intf.fᶠ,a.intf.α,a.intf.n̂,u,u,Δt,a.intf.c̄,a.intf.ρuf,zero(T);perdir=a.flow.perdir)
end
# deal with production term
production!(::Nothing, args...) = nothing
production!(prod::Function, u⁰, u, intf) = @loop intf.f[I] = prod(I,u⁰,u,intf.f) over I in inside(intf.f)

# scalar stress
import WaterLily: S
σ_scalar(I,u;μ=1) = √(0.5*sum(abs2,2μ*S(I,u))) # scalar measure for shear stress
hemolysis(I,u⁰,u,f) = min(f[I]+ifelse(σ_scalar(I,u;μ=1) < 0.36, 0.f0, 1.f0), 1.f0) # make sure we don't overfill

function main(;L=64,Re=250,U=1,mem=Array,T=Float32)
    body = AutoBody((x,t)->√sum(abs2,x-SA_F32[2L,2L])-L÷2)
    TwoPhaseSimulation((8L,4L),(U,0),L;ν=U*L/Re,mem,body,T,InterfaceSDF=(x)->1e8)
end

using CUDA
sim = main(;L=64,mem=CuArray)

@gif for tᵢ in range(0,10,step=0.1)
    println(tᵢ)
    while sim_time(sim) < tᵢ
        sim_step!(sim); passive_advect!(sim;production=hemolysis)
    end
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
    flood(sim.flow.σ,clims=(-5,5),axis=([],false),
          cfill=:blues,legend=false,border=:none); body_plot!(sim)
    vol=@inbounds sum(sim.intf.f[inside(sim.intf.f)])/sim.L^2/4/π
    contour!(Array(sim.intf.f)',levels=[0.5],color=:gold,alpha=0.5,lw=2,title="vol(f)/πR²=$(round(vol,digits=2))")
end