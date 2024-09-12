using WaterLily,StaticArrays,WriteVTK,CUDA

function profile!(sim,profiles;loc=[1,2]) # measure velocity profiles at certain locations
    l = []
    for x ∈ loc.*sim.L
        push!(l,azimuthal_avrg(sim.flow.u[Int(x),:,:,1]))
    end
    push!(profiles,l)
end
Base.hypot(I::CartesianIndex) = √sum(abs2,I.I)
using StatsBase, LinearAlgebra
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
    h(x::T) where T = 4L ≤ x ≤ 5L ? convert(T,√stenosis*L/4*0.5*(1-cos(2π*x/L))) : zero(T)
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
    Simulation((10L,L,L), u_pipe, L; U=one(T), ν=U*L/Re, body, mem, T, exitBC=false)
end

sim = make_sim3D(32;mem=Array)
wr = vtkWriter("ThreeD_Stenosis"; attrib=custom_attrib)
t₀,duration,tstep = sim_time(sim),1,0.1;

# run
@time for tᵢ in range(t₀,t₀+duration;step=tstep)
    sim_step!(sim,tᵢ;remeasure=false,verbose=false)
    write!(wr,sim)
    println("tU/L=",round(tᵢ,digits=4),", Δt=",round(sim.flow.Δt[end],digits=3))
end
close(wr)

# plot velocity profiles
using Plots
stenosis = 0.5
profile = []; xs = [1,2,3,4,4.5,5,6,7,8,9,10];
profile!(sim,profile;loc=xs)
p = plot(title="Instantaneous radial U-velocity profiles",size=(1000,200));
for i ∈ 1:length(profile[1])
    u = 0.5*profile[1][i][2].+xs[i]
    y = (profile[1][i][1].-.5)/sim.L.+0.5
    idx = findall(y .≤ 1) # trim to pipe edge
    plot!(p,vcat(reverse(u[idx]),u[idx]),vcat(y[idx].-0.5,y[idx]),
          color=:black,ls=:dash,label=:none)
end
h(x) = 4 ≤ x ≤ 5 ? √stenosis*0.125*(1-cos(2π*x)) : 0
plot!(p,0:0.01:10,h.(0:0.01:10),color=:black,lw=2,label=:none)
plot!(p,0:0.01:10,1.0.-h.(0:0.01:10),color=:black,lw=2,label=:none)
savefig(p,"radial_velocity_profiles.png")
