using WaterLily,StaticArrays,Plots,CUDA

function make_sim2D(L=32;U=1,Re=500,mem=Array,T=Float32)
    # box sdf, center at c and of half-width w
    function box(x,t;c=SA[0,L],w=SA[0.85f0L,1.5f0L])
        d = abs.(x-c)-w; norm(max.(d,0.f0)) + min(max(d[1],d[2]),0.f0)
    end
    function u_pipe(i,x,t)
        i ≠ 1 && return 0.f0
        x[1] ≥ 2 && return 0.f0 # only apply at inlet
        r = abs(x[2]-L) # move to center of domain
        return r<L/24.f0 ? 2.f0-2.f0.*r^2/(L/24.f0)^2.f0 : 0.f0 # remove radius and add stenosis (and the ghost)
    end
    @inline norm(x) = √sum(abs2,x)
    # triangle sdf
    function triangle(p,t)
        q = SA[0.25f0L,L]
        px = abs(p[1])
        a = p .- q.*clamp(p'*q/(q'*q),0.f0,1.f0)
        b = p .- q.*SA[clamp(px/q[1],0.f0,1.f0),1.f0]
        s = -sign(q[2])
        d = min(SA[a'*a, s*(px*q[2]-p[2]*q[1])],
                SA[b'*b, s*(p[2]-q[2])])
        return -sqrt(d[1]).*sign(d[2])
    end
    # constant rotations params
    α = T(-π/2.f0); R = SA{T}[cos(α) sin(α); -sin(α) cos(α)]
    function map(x,t)
        R*(x-SA[L/20.f0,L])
    end
    # make body, https://www.preprints.org/manuscript/202409.1782/v1
    body = AutoBody(box) - (AutoBody((x,t)->√sum(abs2, x.-SA[L/2.f0,L]) .- L/3.f0) -
                            (AutoBody((x,t)->√sum(abs2, x.-SA[L/2.f0,L]) .- L/4.5f0) -
                             AutoBody(triangle,map))) -
           AutoBody((x,t)->box(x,t;c=SA[0.f0,L],w=SA[2.1f0L,L/24.f0]))

    # make sim
    Simulation((4L,2L), u_pipe, L; U, ν=U*L/Re, body, mem, T, exitBC=true)
end

# perturb
function force!(flow,t)
    @WaterLily.loop flow.f[I,2] += 0.0 over I ∈ CartesianIndices((380:400,440:460))
end

sim = make_sim2D(128;mem=Array)

# quick check that we have what we want
measure_sdf!(sim.flow.σ,sim.body)
flood(sim.flow.σ,clims=(-1,1))

# run
t₀,duration,tstep = sim_time(sim),20,0.2;
@gif for tᵢ in range(t₀,t₀+duration;step=tstep)
    sim_step!(sim,tᵢ;remeasure=false,verbose=false,udf=force!)
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
    @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,abs(sim.flow.σ[I]))
    flood(sim.flow.σ, clims=(0,20), levels=21, axis=([], false), cfill=cgrad(:bone_1, rev=true),
          legend=false,border=:none,size=size(sim.flow.p).-2)
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end