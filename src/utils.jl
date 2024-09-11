using LinearAlgebra: dot, norm
"""
    ellipse(p, ab; itr=5)
    -iter: Number of Newton iterations
"""
function ellipse(xyz, ab)
    # put in symmetry
    x = SA[abs(xyz[3]),√(xyz[1]^2+xyz[2]^2)]

    # find root with Newton solver
    q = ab.*(x-ab);
    w = (q[1]<q[2]) ? π/2 : 0.0;
    for i ∈ 1:5
        u = ab.*SA[ cos(w),sin(w)]
        v = ab.*SA[-sin(w),cos(w)]
        w += dot(x-u,v)/(dot(x-u,u)+dot(v,v));
    end
    # compute final point and distance
    d = norm(x-ab.*SA[cos(w),sin(w)]);
    
    # return signed distance
    return (dot(x./ab,x./ab)>1.0) ? d : -d
end


# make a writer with some attributes
velocity(a::Simulation) = a.flow.u |> Array;
pressure(a::Simulation) = a.flow.p |> Array;
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a.flow)); 
                        a.flow.σ |> Array;)
vorticity(a::Simulation) = (@inside a.flow.σ[I] = 
                            WaterLily.curl(3,I,a.flow.u)*a.L/a.U;
                            a.flow.σ |> Array;)
_vbody(a::Simulation) = a.flow.V |> Array;
mu0(a::Simulation) = a.flow.μ₀ |> Array;
lamda(a::Simulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I, a.flow.u);
                        a.flow.σ |> Array;)
_vbody(a::Simulation) = a.flow.V |> Array;

custom_attrib = Dict(
    "u" => velocity, "p" => pressure, "v" => _vbody,
    "d" => _body, "ω" => vorticity, "λ₂" => lamda,
)# this 