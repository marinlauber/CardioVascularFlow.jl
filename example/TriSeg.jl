using CirculatorySystemModels
using ModelingToolkit

@mtkmodel TriSeg begin
    @components begin
            in = Pin()
            out = Pin()
    end
    @structural_parameters begin
            inP = false
    end
    @variables begin
            V(t) = 2.0
            p(t) = 0.0
    end
    @parameters begin
            V₀
            p₀ = 0.0
            Eₘᵢₙ
            Eₘₐₓ
            n₁
            n₂
            τ
            τ₁
            τ₂
            k
            Eshift = 0.0
    end

    begin
            E = DHelastance(t, Eₘᵢₙ, Eₘₐₓ, n₁, n₂, τ, τ₁, τ₂, Eshift, k)
            DE = DHdelastance(t, Eₘᵢₙ, Eₘₐₓ, n₁, n₂, τ, τ₁, τ₂, Eshift, k)
            p_rel = p₀
    end

    @equations begin
            0 ~ in.p - out.p
            p ~ in.p
            if inP
                    V ~ (p - p_rel) / E + V₀
                    D(p) ~ (in.q + out.q) * E + (p - p_rel) / E * DE
            else
                    p ~ (V - V₀) * E + p_rel
                    D(V) ~ in.q + out.q
            end
    end
end