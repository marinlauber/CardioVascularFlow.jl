using CirculatorySystemModels
using ModelingToolkit
using OrdinaryDiffEq
using Plots
include("TriSeg.jl")

## Start Modelling
@variables t

## Ventricles
@named LV = ShiChamber(V₀=v0_lv, p₀=p0_lv, Eₘᵢₙ=Emin_lv, Eₘₐₓ=Emax_lv, τ=τ, τₑₛ=τes_lv, τₑₚ=τed_lv, Eshift=0.0)
# The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
@named LA = ShiChamber(V₀=v0_la, p₀=p0_la, Eₘᵢₙ=Emin_la, Eₘₐₓ=Emax_la, τ=τ, τₑₛ=τpww_la / 2, τₑₚ=τpww_la, Eshift=τpwb_la)
@named RV = ShiChamber(V₀=v0_rv, p₀=p0_rv, Eₘᵢₙ=Emin_rv, Eₘₐₓ=Emax_rv, τ=τ, τₑₛ=τes_rv, τₑₚ=τed_rv, Eshift=0.0)
# The atrium can be defined either as a ShiChamber with changed timing parameters, or as defined in the paper
# @named RA = ShiChamber(V₀=v0_ra, p₀ = p0_ra, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τₑₛ=τpww_ra/2, τₑₚ =τpww_ra, Eshift=τpwb_ra)
@named RA = ShiAtrium(V₀=v0_ra, p₀=1, Eₘᵢₙ=Emin_ra, Eₘₐₓ=Emax_ra, τ=τ, τpwb=τpwb_ra, τpww=τpww_ra) #, Ev=Inf)

## Valves as simple valves
@named AV = OrificeValve(CQ=CQ_AV)
@named MV = OrificeValve(CQ=CQ_MV)
@named TV = OrificeValve(CQ=CQ_TV)
@named PV = OrificeValve(CQ=CQ_PV)

####### Systemic Loop #######
# Systemic Aortic Sinus ##
@named SAS = CRL(C=Csas, R=Rsas, L=Lsas)
# Systemic Artery ##
@named SAT = CRL(C=Csat, R=Rsat, L=Lsat)
# Systemic Arteriole ##
@named SAR = Resistor(R=Rsar)
# Systemic Capillary ##
@named SCP = Resistor(R=Rscp)
# Systemic Vein ##
@named SVN = CR(R=Rsvn, C=Csvn)

####### Pulmonary Loop #######
# Pulmonary Aortic Sinus ##
@named PAS = CRL(C=Cpas, R=Rpas, L=Lpas)
# Pulmonary Artery ##
@named PAT = CRL(C=Cpat, R=Rpat, L=Lpat)
# Pulmonary Arteriole ##
@named PAR = Resistor(R=Rpar)
# Pulmonary Capillary ##
@named PCP = Resistor(R=Rpcp)
# Pulmonary Vein ##
@named PVN = CR(R=Rpvn, C=Cpvn)

##
circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, SAS.in)
    connect(SAS.out, SAT.in)
    connect(SAT.out, SAR.in)
    connect(SAR.out, SCP.in)
    connect(SCP.out, SVN.in)
    connect(SVN.out, RA.in)
    connect(RA.out, TV.in)
    connect(TV.out, RV.in)
    connect(RV.out, PV.in)
    connect(PV.out, PAS.in)
    connect(PAS.out, PAT.in)
    connect(PAT.out, PAR.in)
    connect(PAR.out, PCP.in)
    connect(PCP.out, PVN.in)
    connect(PVN.out, LA.in)
    connect(LA.out, MV.in)
    connect(MV.out, LV.in)
]

## Compose the whole ODE system
@named _circ_model = ODESystem(circ_eqs, t)
@named circ_model = compose(_circ_model,
    [LV, RV, LA, RA, AV, MV, PV, TV, SAS, SAT, SAR, SCP, SVN, PAS, PAT, PAR, PCP, PVN])

## And simplify it
circ_sys = structural_simplify(circ_model)

## Setup ODE
# Initial Conditions for Shi Valve
# u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas , pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn, 0, 0, 0, 0,0, 0, 0, 0]
# and for OrificeValve --- Commment this next line to use ShiValves
# u0 = [LV_Vt0, RV_Vt0, LA_Vt0, RA_Vt0, pt0sas, qt0sas, pt0sat, qt0sat, pt0svn, pt0pas, qt0pas, pt0pat, qt0pat, pt0pvn]

u0 = [
    LV.V => LV_Vt0
    LV.p => (LV_Vt0 - v0_lv) * Emin_lv + p0_lv
    RV.V => RV_Vt0
    RV.p => (RV_Vt0 - v0_rv) * Emin_rv + p0_rv
    LA.V => LA_Vt0
    RA.V => RA_Vt0
    SAS.C.p => pt0sas
    SAS.C.V => pt0sas * Csas
    SAS.L.q => qt0sas
    SAT.C.p => pt0sat
    SAT.C.V => pt0sat * Csat
    SAT.L.q => qt0sat
    SVN.C.p => pt0svn
    SVN.C.V => pt0svn * Csvn
    PAS.C.p => pt0pas
    PAS.C.V => pt0pas * Cpas
    PAS.L.q => qt0pas
    PAT.C.p => pt0pat
    PAT.C.V => pt0pat * Cpat
    PAT.L.q => qt0pat
    PVN.C.p => pt0pvn
    PVN.C.V => pt0pvn * Cpvn
]

prob = ODEProblem(circ_sys, u0, (0.0, 20.0))
##
@time sol = solve(prob, Tsit5(), reltol=1e-9, abstol=1e-12, saveat=19:0.01:20)

p1 = plot(sol, idxs=[LV.p, RA.p], tspan=(19, 20), xlabel = "Time [s]", ylabel = "Pressure [mmHg]",  hidexaxis = nothing) # Make a line plot
p2 = plot(sol, idxs=[LV.V], tspan=(19, 20),xlabel = "Time [s]", ylabel = "Volume [ml]",  linkaxes = :all)
p3 = plot(sol, idxs=[SAS.L.q,PAS.L.q], tspan=(19, 20),xlabel = "Time [s]", ylabel = "Flow rate [ml/s]", linkaxes = :all)
p4 = plot(sol, idxs=(LV.V, LV.p),xlabel = "Volume [ml]", ylabel = "Pressure [mmHg]", linkaxes = :all)

img = plot(p1, p2, p3, p4; layout=@layout([a b; c d]), legend = true)
img