using WaterLily,StaticArrays,WriteVTK,CUDA

# make a writer with some attributes
vtk_velocity(a::Simulation) = a.flow.u |> Array;
vtk_pressure(a::Simulation) = a.flow.p |> Array;
vtk_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); 
                          a.flow.σ |> Array;)
custom_attrib = Dict("u" => vtk_velocity, "p" => vtk_pressure, "d" =>vtk_body)

# make the sim
function make_pipe(L=32;stenosis=0.5,U=1,Re=2500,mem=Array,T=Float32)
    h(x::T) where T = 3L ≤ x ≤ 4L ? convert(T,√stenosis*L/4*0.5*(1-cos(2π*x/L))) : zero(T)
    function pipe(x,t)
        r = √sum(abs2,SA[x[2],x[3]].-L/2.f0) # move to center of pipe
        L/2.f0 - r - 1.5f0 - h(x[1]) # remove radius and add stenosis (and the ghost)
    end
    # analytical solution laminar pipe flow u/U ~ (1-y^2/R^2) ∀ y ∈ [0,L/2]
    function u_pipe(i,x,t)
        i ≠ 1 && return 0.f0
        r = √sum(abs2,SA[x[2],x[3]].-L/2.f0)
        return r<L/2.f0-1.5f0 ? 2.f0-2.f0.*r^2/(L/2.f0-1.5f0)^2.f0 : 0.f0 # remove radius and add stenosis (and the ghost)
    end
    # pressure gradient required to drive the flow to u~1
    body = AutoBody(pipe)
    Simulation((10L,L,L), u_pipe, L; U=one(T), ν=U*L/Re, body, mem, T, exitBC=true)
end
# make the sim
function make_channel(L=32;stenosis=0.5,U=1,Re=2500,mem=Array,T=Float32)
    h(x::T) where T = 3L ≤ x ≤ 4L ? convert(T,√stenosis*L/4*0.5*(1-cos(2π*x/L))) : zero(T)
    function channel(x,t)
        L/2.f0 - abs(x[2]-L/2) - 1.5f0 - ifelse(x[2]>L/2, 0, h(x[1]))
    end
    # analytical solution laminar channel flow u/U ~ 9(y/L-y^2/L^2) ∀ y ∈ [0,L]
    function u_channel(i,x,t)
        i ≠ 1 && return 0.f0
        r,y = abs(x[2]-L/2), x[2]-1.5f0
        return r<L/2.f0-1.5f0 ? 9.f0*(y/(L.-3.0f0)-(y/(L.-3.0f0))^2.f0) : 0.f0
    end
    # pressure gradient required to drive the flow to u~1
    body = AutoBody(channel)
    Simulation((10L,L,L), u_channel, L; U=one(T), ν=U*L/Re, body, mem, T, perdir=(3,), exitBC=true)
    # g(i,x,t) = i == 1 ? U^2.f0/(L/2.f0)^2.f0 : 0.f0
    # Simulation((10L,L,L), (0,0,0), L; U=one(T), ν=U*L/Re, body, mem, T, g, perdir=(1,3))
end
# fda nozzel sim https://doi.org/10.1007/s11517-020-02188-8
function make_FDA_nozzel(L=32;U=1.f0/9.f0,Re=2500,mem=Array,T=Float32)
    # very strange geometry with anoying ratios
    L₁=L*T(4537/2400); L₂=L*T(10/3); m₁=T(800/4537)
    S(x::T) where T = 0 ≤ x-3L ≤ L₁+L₂ ? convert(T,min((x[1]-3L)*m₁,L/3)) : zero(T)
    function pipe(x,t)
        r = √sum(abs2,SA[x[2],x[3]].-L/2.f0) # move to center of pipe
        L/2.f0 - r - 1.5f0 - S(x[1]) # remove radius and add stenosis (and the ghost)
    end
    # analytical solution laminar pipe flow u/U ~ (1-y^2/R^2) ∀ y ∈ [0,L/2]
    function u_pipe(i,x,t)
        i ≠ 1 && return 0.f0
        r = √sum(abs2,SA[x[2],x[3]].-L/2.f0)
        return r<L/2.f0-1.5f0 ? 2*U*(1.f0-r^2/(L/2.f0-1.5f0)^2.f0) : 0.f0 # remove radius and add stenosis (and the ghost)
    end
    # pressure gradient required to drive the flow to u~1
    body = AutoBody(pipe)
    Simulation((16L,L,L), u_pipe, L; U, ν=U*L/Re, body, mem, T, exitBC=true)
end

# sim = make_channel(48;mem=CuArray)
# sim = make_pipe(48;mem=CuArray)
sim = make_FDA_nozzel(96;mem=CuArray)
wr = vtkWriter("uBC_test"; attrib=custom_attrib)
t₀,duration,tstep = sim_time(sim),50,0.2;

# run
@time for tᵢ in range(t₀,t₀+duration;step=tstep)
    sim_step!(sim,tᵢ;remeasure=false,verbose=false)
    write!(wr,sim)
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr)

using Plots
σ = Array(zeros(size(sim.flow.σ[:,:,1]))); u = Array(sim.flow.u);
J(I) = CartesianIndex(I[1],I[2],size(sim.flow.σ,3)÷2)
@inside σ[I] = ifelse(sim.body.sdf(loc(0,J(I)),0)≥0,√sum(abs2,WaterLily.ω(J(I),u)*sim.L/sim.U),NaN);
flood(σ, clims=(-125,1250/3), axis=([], false), cfill=cgrad(:bone_1, rev=true),
      legend=false,border=:none,size=(10*sim.L,sim.L), dpi=1200)
savefig("fda_nozzel_vorticity.png")

# plot velocity profiles
using Plots,ForcePartition
let
    profile = []; xs = collect(0:1:15);
    for x ∈ xs.*sim.L
        push!(profile,ForcePartition.azimuthal_avrg(u[max(1,Int(x)),:,:,1]))
    end
    p = plot(size=(1000,200), aspect_ratio=:equal);
    for i ∈ 1:length(profile)
        uᵢ = 0.6.*profile[i][2].+xs[i]
        y = (profile[i][1].-0.5)/sim.L.+0.5
        idx = findall(y .≤ 1) # trim to pipe edge
        plot!(p,vcat(reverse(uᵢ[idx]),uᵢ[idx]),vcat(y[idx].-0.5,y[idx]),
              color=:black,label=:none,lw=2)
        plot!(p,[xs[i],xs[i]],[0,1],color=:black,alpha=0.2,lw=0.5,label=:none)
    end
    L₁=4537/2400; L₂=10/3; m₁=800/4537
    S(x::T) where T = 0 ≤ x-3≤ L₁+L₂ ? convert(T,min((x[1]-3)*m₁,1/3)) : zero(T)
    plot!(p,0:0.01:16,S.(0:0.01:16),color=:black,lw=1,label=:none)
    plot!(p,0:0.01:16,1.0.-S.(0:0.01:16),color=:black,lw=1,label=:none)
# savefig(p,"radial_velocity_profiles.png")
end