using CirculatorySystemModels
using ModelingToolkit
using OrdinaryDiffEq
using Plots
# import DisplayAs

# Cycle time in seconds
τ = 0.85

# Double Hill parameters for the ventricle
Eₘᵢₙ = 0.03
Eₘₐₓ = 1.5
n1LV = 1.32
n2LV = 21.9
Tau1fLV = 0.303 * τ
Tau2fLV = 0.508 * τ

# Resistances and Compliances
R_s = 1.11
C_sa = 1.13
C_sv = 11.0

# Aortic valve basic
Zao = 0.033
# Mitral valve basic
Rmv = 0.006

# Inital Pressure (mean cardiac filling pressure)
MCFP = 7.0

# Calculating the additional `k` parameter
t = LinRange(0, τ, 1000)
kLV = 1 / maximum((t ./ Tau1fLV) .^ n1LV ./ (1 .+ (t ./ Tau1fLV) .^ n1LV) .* 1 ./ (1 .+ (t ./ Tau2fLV) .^ n2LV))


# Set up the model elements
@parameters t

# Heart is modelled as a single chamber (we call it `LV` for "Left Ventricle" so the model can be extended later, if required):
@named LV = DHChamber(V₀=0.0, Eₘₐₓ=Eₘₐₓ, Eₘᵢₙ=Eₘᵢₙ, n₁=n1LV, n₂=n2LV, τ=τ, τ₁=Tau1fLV, τ₂=Tau2fLV, k=kLV, Eshift=0.0, inP=true)

# The two valves are simple diodes with a small resistance
# (resistance is needed, since perfect diodes would connect two elastances/compliances, which will lead to unstable oscillations):
@named AV = ResistorDiode(R=Zao)
@named MV = ResistorDiode(R=Rmv)

# The main components of the circuit are 1 resistor `Rs` and two compliances for systemic arteries `Csa`,
# and systemic veins `Csv` (names are arbitrary).
@named Rs = Resistor(R=R_s)
@named Csa = Compliance(C=C_sa, inP=true, has_ep=true, has_variable_ep=true)
@named Csv = Elastance(E=1/C_sv, inP=false)

circ_eqs = [
    connect(LV.out, AV.in)
    connect(AV.out, Csa.in)
    connect(Csa.out, Rs.in)
    connect(Rs.out, Csv.in)
    connect(Csv.out, MV.in)
    connect(MV.out, LV.in)
    Csa.ep.p ~ 0
]

# Add the component equations
@named _circ_model = ODESystem(circ_eqs, t)

@named circ_model = compose(_circ_model,
    [LV, AV, MV, Rs, Csa, Csv])

# Simplify the ODE system
circ_sys = structural_simplify(circ_model)

# initial conditions
u0 = [
    LV.p => MCFP
    LV.V => 10 + (MCFP - 1) / Eₘᵢₙ
    Csa.p => MCFP
    Csa.V => MCFP*C_sa
    Csv.p => MCFP
    Csv.V => MCFP*C_sv
]

# Then we can define the problem:
tspan = (0, 20)
prob = ODEProblem(circ_sys, u0, tspan)

# Simulate
@time sol = solve(prob, Vern7(), reltol=1e-6, abstol=1e-9, saveat=(19:0.01:20))

p1 = plot(sol, idxs=[LV.p,  Csa.in.p], tspan=(19, 20), xlabel = "Time [s]", ylabel = "Pressure [mmHg]",  hidexaxis = nothing) # Make a line plot
p2 = plot(sol, idxs=[LV.V], tspan=(19, 20),xlabel = "Time [s]", ylabel = "Volume [ml]",  linkaxes = :all)
p3 = plot(sol, idxs=[Csa.in.q,Csv.in.q], tspan=(19, 20),xlabel = "Time [s]", ylabel = "Flow rate [ml/s]", linkaxes = :all)
p4 = plot(sol, idxs=(LV.V, LV.p),xlabel = "Volume [ml]", ylabel = "Pressure [mmHg]", linkaxes = :all)

img = plot(p1, p2, p3, p4; layout=@layout([a b; c d]), legend = true)
savefig(img,"single_chamber_model.png")