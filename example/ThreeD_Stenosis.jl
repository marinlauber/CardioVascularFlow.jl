using WaterLily,StaticArrays,WriteVTK,CUDA
include("../src/utils.jl")

# make the sim
function make_sim3D(L=32;stenosis=0.5,U=1,Re=2500,mem=Array,T=Float32)
    h(x::T) where T = 3L ≤ x ≤ 4L ? convert(T,√stenosis*L/4*0.5*(1-cos(2π*x/L))) : zero(T)
    function pipe(x,t)
        r = √sum(abs2,SA[x[2],x[3]].-L/2.f0) # move to center of pipe
        L/2 - r - 1.5f0 - h(x[1]) # remove radius and add stenosis (and the ghost)
    end
    # analytical solution laminar pipe flow u/U ~ (1-y^2/L^2) ∀ y ∈ [0,L/2]
    function u_pipe(i,x,t)
        i ≠ 1 && return 0.f0
        r = √sum(abs2,SA[x[2],x[3]].-L/2.f0)
        ifelse(r<L/2.f0-1.5f0,1.f0-r^2/L^2.f0,0.f0) # remove radius and add stenosis (and the ghost)
    end
    body = AutoBody(pipe)
    Simulation((10L,L,L), u_pipe, L; U=one(T), ν=U*L/Re, body, mem, T, exitBC=false)
end

sim = make_sim3D(64;mem=CuArray)
wr = vtkWriter("ThreeD_Stenosis"; attrib=custom_attrib)
t₀,duration,tstep = sim_time(sim),0.2,0.1;

# run
@time for tᵢ in range(t₀,t₀+duration;step=tstep)
    sim_step!(sim,tᵢ;remeasure=false,verbose=false)
    write!(wr,sim)
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr)

# # plot velocity profiles
# using Plots
# stenosis = 0.5
# profile = []; xs = [1,2,3,4,4.5,5,6,7,8,9,10];
# profile!(sim,profile;loc=xs)
# p = plot(title="Instantaneous radial U-velocity profiles",size=(1000,200));
# for i ∈ 1:length(profile[1])
#     u = 0.5*profile[1][i][2].+xs[i]
#     y = (profile[1][i][1].-.5)/sim.L.+0.5
#     idx = findall(y .≤ 1) # trim to pipe edge
#     plot!(p,vcat(reverse(u[idx]),u[idx]),vcat(y[idx].-0.5,y[idx]),
#           color=:black,ls=:dash,label=:none)
# end
# h(x) = 4 ≤ x ≤ 5 ? √stenosis*0.125*(1-cos(2π*x)) : 0
# plot!(p,0:0.01:10,h.(0:0.01:10),color=:black,lw=2,label=:none)
# plot!(p,0:0.01:10,1.0.-h.(0:0.01:10),color=:black,lw=2,label=:none)
# savefig(p,"radial_velocity_profiles.png")
