using WaterLily
using ForwardDiff
using StaticArrays
import WaterLily: @loop,CI,∂

@inline F(i,j,I,ξ) = inv(∂(j,CI(I,i),ξ)+eps())
# uniform field, no deformation gradient
ξ = ones(10,10,2); apply!((i,x)->x[i], ξ)
σ = ones(10,10)
@inside σ[I] = F(1,1,I,ξ)
@inline FFᵀ(i,j,I,ξ) = inv(∂(j,CI(I,i),ξ))*inv(∂(i,CI(I,j),ξ))
@inline trFFᵀ(i,j,I::CartesianIndex{n},ξ) where n = sum(FFᵀ(i,i,I,ξ) for i ∈ 1:n)
# plane-strain incompressible neo Hookean law
@inline τ(i,j,I,ξ) = (FF(i,j,I,ξ)-ifelse(i==j,inv(3)*(trFFᵀ(i,j,I,ξ)+1),0))

struct ReferenceMapBody{T,F1<:Function,F2<:Function,Vf<:AbstractArray{T}} <: AbstractBody
    sdf :: F1
    map :: F2
    ξ   :: Vf # reference map field
    ξ⁰  :: Vf # initial reference map field
    V   :: Vf # Eulerian displacement field
    function ReferenceMapBody(N::NTuple{D}, sdf, map=(x,t)->x; compose=true, T=Float32, f=Array) where D
        comp(x,t) = compose ? sdf(map(x,t),t) : sdf(x,t)
        Ng = N .+ 2; Nd = (Ng..., D)
        ξ = Array{T}(undef, Nd...) |> f; apply!((i,x)->x[i], ξ)  # ξ(x,0) = ξ⁰(x) = x
        ξ⁰= Array{T}(undef, Nd...) |> f; apply!((i,x)->x[i], ξ⁰)
        V = Array{T}(undef, Nd...) |> f; V .= 0.0
        new{T,typeof(comp),typeof(map),typeof(ξ)}(comp, map, ξ, ξ⁰, V)
    end
end

WaterLily.sdf(body::ReferenceMapBody,x,t;kwargs...) = body.sdf(x,t)

function WaterLily.measure(body::ReferenceMapBody,x,t;fastd²=Inf)
    # eval d=f(x,t), and n̂ = ∇f
    d = body.sdf(x,t)
    d^2>fastd² && return (d,zero(x),zero(x)) # skip n,V
    n = ForwardDiff.gradient(x->body.sdf(x,t), x)
    any(isnan.(n)) && return (d,zero(x),zero(x))

    # correct general implicit fnc f(x₀)=0 to be a pseudo-sdf
    #   f(x) = f(x₀)+d|∇f|+O(d²) ∴  d ≈ f(x)/|∇f|
    m = √sum(abs2,n); d /= m; n /= m

    # The velocity depends on the material change of ξ=m(x,t):
    #   Dm/Dt=0 → ṁ + (dm/dx)ẋ = 0 ∴  ẋ =-(dm/dx)\ṁ
    J = ForwardDiff.jacobian(x->body.map(x,t), x)
    dot = ForwardDiff.derivative(t->body.map(x,t), t)
    return (d,n,-J\dot)
end

function body_update!(V,V⁰,r,ξ,Φ,dt;G=32)
    r .= 0.
    N,n = size_u(V)
    for i ∈ 1:n, j ∈ 1:n
        # compute τ (∇ξ⁻¹∇ξ⁻ᵀ) at cell center and add ∇⋅τ (∇⋅∇ξ⁻¹) to the face
        @loop (Φ[I] = G*τ(i,j,I,ξ); r[I,i] += Φ[I]) over I ∈ inside_u(N,j);
        # add the neighboring ∇⋅τ to the face (this finalizes the computation of ∇⋅τ)
        @loop r[I-δ(j,I),i] -= Φ[I] over I ∈ inside_u(N,j)
    end
    # update reference map velocity
    @loop V[Ii] = V⁰[Ii]+dt*r[Ii] over Ii in CartesianIndices(r)
end

function reference_map!(r,u,Φ,ξ⁰,ξ,dt)
    r .= 0.
    N,n = size_u(u)
    #update the reference map ξ
    for i ∈ 1:n, j ∈ 1:n
        # inner cells convection
        @loop (Φ[I] = ϕu(j,CI(I,i),u,ϕ(i,CI(I,j),ξ));
               r[I,i] += Φ[I]) over I ∈ inside_u(N,j)
        @loop r[I-δ(j,I),i] -= Φ[I] over I ∈ inside_u(N,j)
    end
    # update reference map ξ
    @loop ξ[Ii] = ξ⁰[Ii]+dt*r[Ii] over Ii in CartesianIndices(r)
end

import WaterLily: BCTuple,BDIM!,project!,BC!,exitBC!,conv_diff!,accelerate!,CFL,scale_u!
@fastmath function mom_step!(a::Flow{N},body::ReferenceMapBody,b::AbstractPoisson) where N
    a.u⁰ .= a.u; scale_u!(a,0)
    body.ξ⁰ .= body.ξ; scale_ξ!(body,0)
    # predictor u → u'
    U = BCTuple(a.U,@view(a.Δt[1:end-1]),N)
    conv_diff!(a.f,a.u⁰,a.σ,ν=a.ν,perdir=a.perdir)
    body_update!(a.V,body.V,a.f,body.ξ⁰,a.σ,a.Δt[end])
    reference_map!(a.f,a.u⁰,a.σ,body.ξ⁰,body.ξ)
    accelerate!(a.f,@view(a.Δt[1:end-1]),a.g,a.U)
    BDIM!(a); BC!(a.u,U,a.exitBC,a.perdir)
    a.exitBC && exitBC!(a.u,a.u⁰,U,a.Δt[end]) # convective exit
    project!(a,b); BC!(a.u,U,a.exitBC,a.perdir)
    # corrector u → u¹
    U = BCTuple(a.U,a.Δt,N)
    conv_diff!(a.f,a.u,a.σ,ν=a.ν,perdir=a.perdir)
    body_update!(a.V,body.V,a.f,body.ξ,a.σ,a.Δt[end])
    reference_map!(a.f,a.u,a.σ,body.ξ⁰,body.ξ)
    accelerate!(a.f,a.Δt,a.g,a.U)
    scale_ξ!(body,0.5)
    BDIM!(a); scale_u!(a,0.5); BC!(a.u,U,a.exitBC,a.perdir)
    project!(a,b,0.5); BC!(a.u,U,a.exitBC,a.perdir)
    push!(a.Δt,CFL(a))
end
scale_ξ!(a,scale) = @loop a.ξ[Ii] *= scale over Ii ∈ inside_u(front(size(a.ξ)))

function circle(n,m;Re=250,U=1)
    radius, center = m/8, m/2
    body = ReferenceMapBody((n,m), (x,t)->√sum(abs2, x .- center) - radius)
    Simulation((n,m), (U,0), radius; ν=U*radius/Re, body)
end

using Plots
# Initialize and run
sim = circle(3*2^6,2^7)
sim_gif!(sim,duration=10,clims=(-5,5),plotbody=true)

# compute deformation gradient, should be identity tensor
# @SArray [F(i,j,CartesianIndex(10,10),sim.body.ξ) for i ∈ 1:2, j ∈ 1:2]
# @inside sim.flow.σ[I] = F(2,1,I,sim.body.ξ)
# flood(sim.flow.σ,shift=(-2,-1.5),cfill=:seismic)