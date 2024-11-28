using WaterLily
using StaticArrays
include("/home/marin/Workspace/Tutorials-WaterLily/src/TwoD_plots.jl")

# velocity magnitude
mag(I,u) = √sum(ntuple(i->0.25*(u[I,i]+u[I+δ(i,I),i])^2,length(I)))

struct PassiveScalar{T,Sf<:AbstractArray{T}}
    φ⁰:: Sf
    φ :: Sf
    function PassiveScalar(σ)
        new{eltype(σ),typeof(σ)}(copy(σ),copy(σ))
    end
end
import WaterLily: @loop,size_u,ϕ,slice
function update!(s::PassiveScalar,a::Flow;Γ=0.0)
    s.φ⁰ .= s.φ
    s.φ .= 0
    N,n = size_u(a.u)
    for j ∈ 1:n
        # inner cells, cell centered scalar field 
        @loop (a.σ[I] = -a.u[I+δ(j,I),j]*ϕ(j,I+δ(j,I),s.φ⁰);
               s.φ[I] += a.σ[I]) over I ∈ inside(s.φ)
        @loop s.φ[I+δ(j,I)] -= a.σ[I] over I ∈ inside(s.φ)
    end
    # avect
    dt = a.Δt[end]
    @loop s.φ[I] = s.φ⁰[I] + dt*s.φ[I] over I ∈ inside(s.φ)
end

# new BC function
function BC_lid!(a)
    N,n = WaterLily.size_u(a)
    for j ∈ 1:n, i ∈ 1:n
        if i==1 && j==2 # lid, Dirichlet cannot be imposed, must interpolate u[i,j+1,1]+u[i,j]/2 = uBC
            @loop a[I,i] = 2.00-a[I-δ(j,I),i] over I ∈ slice(N,N[j],j)
            @loop a[I,i] = -a[I+δ(j,I),i] over I ∈ slice(N,1,j)
        elseif i==j # Normal direction, homoheneous Dirichlet
            @loop a[I,i] = 0.0 over I ∈ slice(N,1,j)
            @loop a[I,i] = 0.0 over I ∈ slice(N,N[j],j)
        else  # Tangential directions, interpolate ghost cell to homogeneous Dirichlet
            @loop a[I,i] = -a[I+δ(j,I),i] over I ∈ slice(N,1,j)
            @loop a[I,i] = -a[I-δ(j,I),i] over I ∈ slice(N,N[j],j)
        end
    end
end
@fastmath function WaterLily.mom_step!(a::Flow{N},b::AbstractPoisson) where N
    a.u⁰ .= a.u; WaterLily.scale_u!(a,0)
    # predictor u → u'
    WaterLily.conv_diff!(a.f,a.u⁰,a.σ,ν=a.ν,perdir=a.perdir)
    WaterLily.BDIM!(a); BC_lid!(a.u)
    WaterLily.project!(a,b); BC_lid!(a.u)
    # corrector u → u¹
    WaterLily.conv_diff!(a.f,a.u,a.σ,ν=a.ν,perdir=a.perdir)
    WaterLily.BDIM!(a); WaterLily.scale_u!(a,0.5); BC_lid!(a.u)
    WaterLily.project!(a,b,0.5); BC_lid!(a.u)
    push!(a.Δt,WaterLily.CFL(a))
end

# Initialize simulation
L = 2^7
U = 1.0
Re = 100
# using CUDA
sim = Simulation((L,L),(0.0,0.0),L;U=U,ν=U*L/Re,mem=Array)
scalar = PassiveScalar(sim.flow.σ)
apply!(x->ifelse(√sum(abs2,x.-SA[L/4,L/4])-L÷8<0,1.,0.),scalar.φ)

# get start time
t₀ = round(sim_time(sim))
duration = 10; step = 0.1

anim = @animate for tᵢ in range(t₀,t₀+duration;step)

    # update until time tᵢ in the background
    # sim_step!(sim,tᵢ)
    while sim_time(sim) < tᵢ
        mom_step!(sim.flow,sim.pois)
        update!(scalar,sim.flow)
    end

    # flood plot
    # @inside sim.flow.σ[I] = mag(I,sim.flow.u)
    # flood(sim.flow.σ|>Array; shift=(-0.5,-0.5),clims=(0,1))
    flood(scalar.φ|>Array; shift=(-0.5,-0.5),clims=(0,1))

    # print time step
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
gif(anim,"lid_cavity.gif",fps=15)
