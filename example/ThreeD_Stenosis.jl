using WaterLily,StaticArrays,WriteVTK,CUDA

# make a writer with some attributes
vtk_velocity(a::Simulation) = a.flow.u |> Array;
vtk_pressure(a::Simulation) = a.flow.p |> Array;
vtk_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); 
                        a.flow.σ |> Array;)
vtk_vorticity(a::Simulation) = (@inside a.flow.σ[I] = 
                            WaterLily.curl(3,I,a.flow.u)*a.L/a.U;
                            a.flow.σ |> Array;)
vtk_vbody(a::Simulation) = a.flow.V |> Array;
vtk_mu0(a::Simulation) = a.flow.μ₀ |> Array;
vtk_lamda(a::Simulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I, a.flow.u);
                        a.flow.σ |> Array;)

custom_attrib = Dict(
    "u" => vtk_velocity, "p" => vtk_pressure, "v" =>vtk_vbody,
    "d" => vtk_body, "ω" => vtk_vorticity, "λ₂" => vtk_lamda,
)# this

# make the sim, following Varghese et al. (2007 )J. Fluid Mech. vol. 582, pp. 253–280 doi:10.1017/S0022112007005848
function make_sim3D(L=32;stenosis=0.25,U=1,Re=2500,mem=Array,T=Float32)
    h(x::T) where T = 2L ≤ x ≤ 4L ? convert(T,L/2*√stenosis*(1-cos(2π*x/2L))/2) : zero(T)
    function pipe(x,t)
        r = √sum(abs2,SA[x[2],x[3]].-L/2.f0) # move to center of pipe
        L/2 - r - 1.5f0 - h(x[1]) # remove radius and add stenosis (and the ghost)
    end
    # analytical solution laminar pipe flow u/U ~ (1-y^2/R^2) ∀ y ∈ [0,L/2]
    function u_pipe(i,x,t)
        i ≠ 1 && return 0.f0
        r = √sum(abs2,SA[x[2],x[3]].-L/2.f0)
        return r<L/2.f0-1.5f0 ? 2.f0-2.f0.*r^2/(L/2.f0-1.5f0)^2.f0 : 0.f0 # remove radius and add stenosis (and the ghost)
    end
    # pressure gradient required to drive the flow to u~1
    body = AutoBody(pipe)
    Simulation((10L,L,L), u_pipe, L; U=one(T), ν=U*L/Re, body, mem, T, exitBC=false)
end

sim = make_sim3D(96;mem=CuArray)
wr = vtkWriter("ThreeD_Stenosis"; attrib=custom_attrib)
t₀,duration,tstep = sim_time(sim),50,0.1;

# run
@time for tᵢ in range(t₀,t₀+duration;step=tstep)
    sim_step!(sim,tᵢ;remeasure=false,verbose=false)
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
write!(wr,sim)
close(wr)

# σ = Array(zeros(size(sim.flow.σ[:,:,1])));
# J(I) = CartesianIndex(I[1],I[2],size(sim.flow.σ,3)÷2)
# @inside σ[I] = ifelse(sim.body.sdf(loc(0,J(I)),0)≥0,√sum(abs2,WaterLily.ω(J(I),u)*sim.L/sim.U),NaN);
# flood(σ, clims=(-125,1250), axis=([], false), cfill=cgrad(:bone_1, rev=true),
#       legend=false,border=:none,size=(10*sim.L,sim.L), dpi=1200)
# savefig("pipe_vorticity.png")

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
