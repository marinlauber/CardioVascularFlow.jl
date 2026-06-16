using WaterLily,StaticArrays,Plots
import WaterLily: @loop
import InterfaceAdvection: BCf!

# Simulation carrying a persistent passive scalar `f` (e.g. an accumulating
# hemolysis dose) transported by the flow, plus a scratch buffer `f⁰` that holds
# fⁿ during the transport step. Both are scalar fields including the ghost layer.
mutable struct ScalarSimulation <: AbstractSimulation
    sim :: Simulation
    f  :: AbstractArray   # passive scalar, persists across steps
    f⁰ :: AbstractArray   # scratch buffer (holds fⁿ during transport)
    function ScalarSimulation(dims, args...; f0=0, mem=Array, kwargs...)
        sim = Simulation(dims, args...; mem, kwargs...)
        T = eltype(sim.flow.u)
        f  = fill(T(f0), dims.+2) |> mem   # scalar field incl. ghost layer
        f⁰ = zeros(T, dims.+2)    |> mem
        new(sim, f, f⁰)
    end
end
# overload properties
Base.getproperty(f::ScalarSimulation, s::Symbol) = s in propertynames(f) ? getfield(f, s) : getfield(f.sim, s)
Base.setproperty!(f::ScalarSimulation, s::Symbol, x) = s in propertynames(f) ? setfield!(f,s,x) : setproperty!(f.sim,s,x)

function passive_advect!(a::ScalarSimulation;production=nothing)
    u⁰,u,Δt = a.flow.u⁰,a.flow.u,a.flow.Δt[end-1]
    # add the production (source) term, then transport the scalar with a bounded
    # (TVD) flux-limited scheme so `f` is free to grow outside [0,1] (e.g. an
    # accumulating hemolysis dose). `f⁰` holds fⁿ during the transport.
    production!(production,u⁰,u,a.f,Δt)
    scalar_transport!(a.f,a.f⁰,u,u⁰,Δt;perdir=a.flow.perdir)
end
# deal with production term (Δt is threaded through so the source can be a rate)
production!(::Nothing, args...) = nothing
production!(prod::Function, u⁰, u, f, Δt) = @loop f[I] = prod(I,u⁰,u,f,Δt) over I in inside(f)

"""
    scalar_transport!(f,buf,u,u⁰,Δt;perdir=())

Transport the passive scalar `f` over `Δt`, solving ∂f/∂t + u·∇f = 0 with a
conservative, van-Leer–limited (TVD) flux scheme using the time-centred face
velocity ū=(u+u⁰)/2. `buf` is a same-size scratch array that holds fⁿ. Unlike
the VOF/PLIC scheme this is bounded but range-unrestricted, so `f` may exceed
[0,1]. The f·(∇·ū) term keeps a uniform field exactly preserved.
"""
function scalar_transport!(f,buf,u,u⁰,Δt;perdir=())
    buf .= f                                                 # keep fⁿ
    @loop f[I] = f[I] - Δt*scalarDiv(I,buf,u,u⁰) over I ∈ inside(f)
    BCf!(f;perdir)
end
@inline function scalarDiv(I,f,u,u⁰)
    s = zero(eltype(f)); d = zero(eltype(f))
    for a in 1:ndims(f)
        Ir = I+δ(a,I)
        uL = 0.5f0*(u[I,a] +u⁰[I,a])                         # lower-face velocity
        uR = 0.5f0*(u[Ir,a]+u⁰[Ir,a])                        # upper-face velocity
        s += scalarFlux(a,Ir,f,uR) - scalarFlux(a,I,f,uL)    # ∇·(ūf)
        d += uR - uL                                         # ∇·ū
    end
    @inbounds s - f[I]*d                                     # → advective form
end
# van Leer limiter and the TVD face value × face velocity through face `If`
@inline ψvanLeer(r) = (r+abs(r))/(1+abs(r))
@inline function scalarFlux(a,If,f,uf)
    up  = uf > 0
    Iu  = up ? If-δ(a,If) : If                               # upwind cell
    Id  = up ? If          : If-δ(a,If)                      # downwind cell
    Iuu = up ? Iu-δ(a,If)  : Iu+δ(a,If)                      # far-upwind cell
    fu = @inbounds f[Iu]; fd = @inbounds f[Id]
    if 1 ≤ Iuu[a] ≤ size(f,a)                                # 2nd-order interior
        fuu = @inbounds f[Iuu]
        Δu  = fu - fuu
        r   = Δu == 0 ? zero(Δu) : (fd-fu)/Δu
        ff  = fu + 0.5f0*ψvanLeer(r)*Δu
    else
        ff  = fu                                             # 1st-order at boundary
    end
    uf*ff
end

import WaterLily: S
σ_scalar(I,u;μ=1) = √(0.5f0*sum(abs2,2μ*S(I,u))) # nondim strain-rate magnitude γ̇·(L/U) (μ=1)

# Dose-based hemolysis (power-law, Giersiepen/Heuser) in the Eulerian linearized
# form (Garon & Farinas, Grigioni): transport a linearized damage scalar `f` with
# source rate (C·τ^α)^(1/β) integrated over the step; HI = f^β at the end.
αhl, βhl, Cphys = 2.416f0, 0.785f0, 3.62f-7   # Pa^-α·s^-β  (Giersiepen, HI as fraction)
# physical reference scales
ρref = 1060f0    # blood density   [kg/m³]
Uref = 1f0       # velocity scale  [m/s]
Lref = 2f-2      # length scale    [m]
Cstar= Cphys*(ρref*Uref^2)^αhl*(Lref/Uref)^βhl

# scale hemolysis from sim dims
function hemolysis(sim;Cstar,αhl,βhl)
    T = eltype(sim.flow.u)
    ν,U,L = T(sim.flow.ν), T(sim.U), T(sim.L)
    (I,u⁰,u,f,Δt) -> f[I] + (Cstar*(ν*σ_scalar(I,u)/U^2)^αhl)^(1f0/βhl)*(U/L)*Δt
end

function main(;L=64,Re=250,U=1,mem=Array,T=Float32)
    body = AutoBody((x,t)->√sum(abs2,x-SA_F32[2L,2L])-L÷2)
    ScalarSimulation((8L,4L),(U,0),L;ν=U*L/Re,mem,body,T)
end

using CUDA
sim = main(;L=64,mem=CuArray)
hemo = hemolysis(sim;Cstar,αhl,βhl)   # scaled based in sim

@gif for tᵢ in range(0,10,step=0.1)
    println(tᵢ)
    while sim_time(sim) < tᵢ
        sim_step!(sim); passive_advect!(sim;production=hemo)
    end
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
    flood(sim.flow.σ,clims=(-5,5),axis=([],false),
          cfill=:blues,legend=false,border=:none); body_plot!(sim)
    # vol=@inbounds sum(sim.f[inside(sim.f)])/sim.L^2/4/π
    contour!((Array(sim.f).^βhl)',levels=[2e-4],color=:gold,alpha=0.5,lw=2,
             title="Contour of HI<0.02%")
end