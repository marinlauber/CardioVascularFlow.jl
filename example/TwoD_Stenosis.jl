using WaterLily,StaticArrays,WriteVTK,Plots
include("/home/marin/Workspace/Tutorials-WaterLily/src/TwoD_plots.jl")

function make_sim2D(L=32;stenosis=0.5,U=1,Re=500,mem=Array,T=Float32)
    h(x::T) where T = 5L ≤ x ≤ 6L ? convert(T,√stenosis*L/4*0.5*(1-cos(2π*x/L))) : zero(T)
    function sdf(x,t)
        r = abs(x[2]-L/2) # move to center of pipe
        L/2 - r - 1.5f0 - h(x[1]) # remove radius and add stenosis (and the ghost)
    end
    body = AutoBody(sdf)
    function u_pipe(i,x,t)
        i ≠ 1 && return zero(T)
        r = x[2]-L/2 # move to center of pipe
        convert(T,ifelse(r<L/2-1.5f0,1-r^2/L^2,0)) # remove radius and add stenosis (and the ghost)
    end
    Simulation((20L,L), u_pipe, L; U, ν=U*L/Re, body, mem, T, exitBC=false)
end

function profile!(sim,profiles) # measure velocity profiles at certain locations
    l = []
    for x ∈ [1,2,3,4,5,5.5,6,7,8,9,10].*sim.L
        push!(l,sim.flow.u[Int(x),:,1])
    end
    push!(profiles,l)
end

sim = make_sim2D(128) #;mem=CuArray)
t₀,duration,tstep = sim_time(sim),10,0.1;

# run
@gif for tᵢ in range(t₀,t₀+duration;step=tstep)
    sim_step!(sim,tᵢ;remeasure=false,verbose=false)
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
    @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,sim.flow.σ[I])
    flood(sim.flow.σ,shift=(-2,-1.5),clims=(-32,32), axis=([], false),
          cfill=:seismic,legend=false,border=:none,size=(1000,200))
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end

# # plot velocity profiles
# stenosis = 0.5
# profile = []; xs = [1,2,3,4,5,5.5,6,7,8,9,10];
# profile!(sim,profile)
# p = plot(size=(1000,200));
# for i ∈ 1:length(profile[1])
#     plot!(p,0.5*profile[1][i][3:end-2].+xs[i],collect(0:1/(length(profile[1][i])-5):1),
#           color=:black,ls=:dash,label=:none)
# end
# h(x) = 5 ≤ x ≤ 6 ? √stenosis*0.125*(1-cos(2π*x)) : 0
# plot!(p,0:0.01:10,h.(0:0.01:10),color=:black,lw=2,label=:none)
# plot!(p,0:0.01:10,1.0.-h.(0:0.01:10),color=:black,lw=2,label=:none)
# savefig(p,"velocity_profiles.png")