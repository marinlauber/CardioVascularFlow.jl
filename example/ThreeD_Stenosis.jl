using WaterLily,StaticArrays,WriteVTK,CUDA

function profile!(sim,profiles) # measure velocity profiles at certain locations
    l = []
    for x ∈ [1,2,3,4,5,5.5,6,7,8,9,10].*sim.L
        push!(l,azimuthal_avrg(x,sim.flow.u))
    end
    push!(profiles,l)
end

function azimuthal_avrg(i,u)
    k = size(u,3)
    return u[Int(i),:,Int(k÷2),1]
end

# make a writer with some attributes
velocity(a::Simulation) = a.flow.u |> Array;
pressure(a::Simulation) = a.flow.p |> Array;
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body); 
                        a.flow.σ |> Array;)
vorticity(a::Simulation) = (@inside a.flow.σ[I] = 
                            WaterLily.curl(3,I,a.flow.u)*a.L/a.U;
                            a.flow.σ |> Array;)
_vbody(a::Simulation) = a.flow.V |> Array;
mu0(a::Simulation) = a.flow.μ₀ |> Array;
lamda(a::Simulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I, a.flow.u);
                        a.flow.σ |> Array;)

custom_attrib = Dict(
    "u" => velocity, "p" => pressure,
    "d" => _body, "ω" => vorticity, "λ₂" => lamda,
)# this maps what to write to the variable name in the .pvd file


# make the sim
function make_sim3D(L=32;stenosis=0.5,U=1,Re=2500,mem=Array,T=Float32)
    h(x::T) where T = 5L ≤ x ≤ 6L ? convert(T,√stenosis*L/4*0.5*(1-cos(2π*x/L))) : zero(T)
    function pipe(x,t)
        r = √sum(abs2,SA[x[2],x[3]].-L/2) # move to center of pipe
        L/2 - r - 1.5f0 - h(x[1]) # remove radius and add stenosis (and the ghost)
    end
    # analytical solution laminar pipe flow u/U ~ 2*(1-y^2/L^2) ∀ y ∈ [0,L/2]
    function u_pipe(i,x,t)
        i ≠ 1 && return zero(T)
        r = √sum(abs2,SA[x[2],x[3]].-L/2) # move to center of pipe
        convert(T,ifelse(r<L/2-1.5f0,1-r^2/L^2,0)) # remove radius and add stenosis (and the ghost)
    end
    body = AutoBody(pipe)
    Simulation((20L,L,L), u_pipe, L; U=one(T), ν=U*L/Re, body, mem, T, exitBC=false)
end

sim = make_sim3D(64;mem=Array)
wr = vtkWriter("ThreeD_Stenosis"; attrib=custom_attrib)
t₀,duration,tstep = sim_time(sim),100,0.1;

# run
@time for tᵢ in range(t₀,t₀+duration;step=tstep)
    sim_step!(sim,tᵢ;remeasure=false,verbose=false)
    write!(wr,sim)
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr)

# plot velocity profiles
# using Plots
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
# p