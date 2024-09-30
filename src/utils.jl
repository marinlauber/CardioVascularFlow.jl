using LinearAlgebra: dot, norm
using StatsBase
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

function profile!(sim,profiles;loc=[1,2]) # measure velocity profiles at certain locations
    l = []
    for x ∈ loc.*sim.L
        push!(l,azimuthal_avrg(sim.flow.u[Int(x),:,:,1]))
    end
    push!(profiles,l)
end
Base.hypot(I::CartesianIndex) = √sum(abs2,I.I)
function azimuthal_avrg(data; center=nothing, binsize=1.0)
    """
    Calculate the azimuthally averaged radial profile.
    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fractional pixels).
    binsize - size of the averaging bin.  Can lead to strange results if
        non-binsize factors are used to specify the center and the binsize is
        too large
    """
    # Calculate the indices from the image
    CIs = CartesianIndices(data)
    isnothing(center) && (center = (maximum(CIs)-minimum(CIs)).I.÷2)
    
    # radial distance from the center, make it a vector
    r = weights(hypot.(collect(CIs .- CartesianIndex(center)))).values
    
    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)  
    nbins = Int(round(maximum(r) / binsize)+1)
    maxbin = nbins * binsize
    bins = range(0,maxbin;length=nbins+1)
    # but we're probably more interested in the bin centers than their left or right sides...
    bin_centers = (bins[1:end-1].+bins[2:end])/2.0
    r_weights = fit(Histogram, r, bins, closed=:left).weights

    # compute the azimuthal average
    radial_prof = fit(Histogram, r, weights(data), bins, closed=:left).weights ./ r_weights
    
    return bin_centers, radial_prof
end