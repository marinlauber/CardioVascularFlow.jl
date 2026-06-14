using WaterLily,StaticArrays,Plots,OrdinaryDiffEq


function make_sim2D(L=32;stenosis=0.5,U=1,Re=500,mem=Array,T=Float32)

    # Cardiac parameters
    Emax = 2                    #mmHg/ml; slope of the ESPVR
    V0 = 20                     #ml; intercept with volume axis of the ESPVR
    Pfilling = 5                #mmHg; venous filling pressure
    Pout = 5                    #mmHg: 
    EDV = 120                   #ml; end-diastolic volume. We will use EDV with Pvenous to calculate Emin
    Emin = Pfilling/(EDV-V0)    #mmHg/ml
    HR = 60                     #heart rate in beats/min

    # Valve resistances
    Rmv_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
    Rmv_bwd = 1e10              #mmHg/ml/s; leak resistance
    Rao_fwd = 0.002             #mmHg/ml/s; resistance in forward flow direction
    Rao_bwd = 1e10              #mmHg/ml/s; leak resistance

    # Arterial model parameters
    R_WK2 = 1                   #mmHg/ml/s
    C_WK2 = 2                   #ml/mmHg

    # lump together
    params = [Pfilling, Rmv_fwd, Rmv_bwd, Rao_fwd,
              Rao_bwd, R_WK2, C_WK2, Pout]

    # Double Hill function inspired by Stergiopulos et al. (DOI:10.1152/ajpheart.1996.270.6.H2050)
    function Elastance(t;Emin=0.05,Emax=2,a₁=0.303,a₂=0.508,n₁=1.32,n₂=21.9,α=1.672)
        (Emax-Emin) * α * (t%1/a₁)^n₁ / (1+(t%1/a₁)^n₁) * inv(1+(t%1/(a₂))^n₂)  + Emin
    end
    @inline computePLV(t,V;Emin=0.05,Emax=2,V0=20) = Elastance(t;Emin,Emax) * (V-V0)


    function Windkessel!(du,u,p,t)
        # unpack
        (VLV,Pao,_,_) = u
        (Pfill,Rmv_fwd,Rmv_bwd,Rao_fwd,Rao_bwd,R,C,Pout)  = p

        # first calculate PLV from elastance and VLV 
        PLV = computePLV(t,VLV)

        # calculate Qmv, Pfilling>PLV; forward transmitral flow, PLV>Pfilling - backward transmitral flow
        Qmv = Pfill ≥ PLV ? (Pfill-PLV)/Rmv_fwd : (PLV-Pfill)/Rmv_bwd
        u[3] = Qmv # store

        #calculate Qao, PLV>Pao; forward aortic flow, PLV>Pao; backward aortic flow
        Qao = PLV ≥ Pao ? (PLV-Pao)/Rao_fwd : (Pao-PLV)/Rao_bwd
        u[4] = Qao # store

        # rates
        du[1] = Qmv - Qao                 #dVLV/dt=Qmv-Qao
        du[2] = Qao/C - (Pao-Pout)/(R*C)  #dPao/dt=Qao/C-Pao/RC
        du[3] = 0                         # no rate of change
        du[4] = 0                         # no rate of change
    end

    #Setup
    u₀ = [EDV, 60, 0, 0] # initial conditions
    tspan = (0.0, 20.0)

    #Pass to solver
    prob = ODEProblem(Windkessel!, u₀, tspan, params)
    integrator = init(prob, Tsit5(), dtmax=1e-3, reltol=1e-6, abstol=1e-9,
                             save_everystep=false)

    h(x::T) where T = 5L ≤ x ≤ 6L ? convert(T,√stenosis*L/4*0.5*(1-cos(2π*x/L))) : zero(T)
    function sdf(x,t)
        r = abs(x[2]-L/2) # move to center of pipe
        L/2 - r - 3/2 - h(x[1]) # remove radius and add stenosis (and the ghost)
    end
    function u_pipe(i,x,t)
        i ≠ 1 && return zero(T)
        r = abs(x[2]-L/2) # move to center of pipe
        dt = t - integrator.t # how much are we lacking behind
        # @show t, integrator.t, dt
        convert(T,ifelse(r<L/2-3/2,2U*(1-r^2/(L/2-3/2)^2),0)) # remove radius and add stenosis (and the ghost)
    end
    Simulation((20L,L), u_pipe, L; U, ν=U*L/Re, body=AutoBody(sdf), mem, T, exitBC=true)
end

function profile!(sim,profiles) # measure velocity profiles at certain locations
    l = []
    for x ∈ [1,2,3,4,5,5.5,6,7,8,9,10].*sim.L
        push!(l,sim.flow.u[Int(x),:,1])
    end
    push!(profiles,l)
end

sim = make_sim2D(64) #;mem=CuArray)
t₀,duration,tstep = sim_time(sim),1,0.1;

# run
@gif for tᵢ in range(t₀,t₀+duration;step=tstep)
    sim_step!(sim,tᵢ;remeasure=false,verbose=false)
    @inside sim.flow.σ[I] = WaterLily.curl(3,I,sim.flow.u)*sim.L/sim.U
    @inside sim.flow.σ[I] = ifelse(abs(sim.flow.σ[I])<0.001,0.0,abs(sim.flow.σ[I]))
    flood(sim.flow.σ,clims=(0,32), axis=([], false), cfill=cgrad(:bone_1, rev=true),
          legend=false,border=:none,size=(1000,200))
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end