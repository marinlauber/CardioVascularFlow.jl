using StaticArrys
using LinearAlgebra

struct Valve

end

struct Periphery
    
end

struct CircAdapt{T}
    ρ :: T
    ts :: T
    function CircAdapt(;T=Float32) where T
        ρ = 1050.0
        ts = 0.585
        # Tubes: aorta (AO), arteria pulmonalis (AP), venae cavae (VC), and venae pulmonales (VP)
        Atwall = SA{T}[274 242 58 85]
        lt = SA{T}[500 200 400 200]
        kt = SA{T}[5.0 8.0 10.0 10.0]
        # Chambers: left (LV) and right (RV) ventricle; left (LA) and right atrium (RA)
        Vc = SA{T}[57.0 75.3 44.2 54.4]
        hc = SA{T}[15.0 4.0 2.0 2.0]
        Δtact = SA{T}[0.1 0.1 0.02 0.0]
        # Valves: aortic (AV), pulmonary (PV), mitral (MV), and tricuspid (TV) valve; pulmonary (PO) and systemic (SO) outlet
        Aopen = SA{T}[400 400 500 500 400 400]
        Aleak = SA{T}[0.0 0.0 0.0 0.0 400 400] 
        # Periphery: systemic (sys) and pulmonary (pulm) circulation
        Δpref = SA{T}[1.5 10.0]
        qref = SA{T}[85 85]
        rpy = SA{T}[1 2]
        new{T}(ρ, ts)
    end
end

function step!(::CircAdapt{T}) where T
    nothing
end