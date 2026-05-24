using WaterLily
using StaticArrays
using WriteVTK
include("../src/utils.jl")

function leaflet(L=2^5;Re=5e2,mem=Array,U=1,AR=2,T=Float32)
    # Define simulation size, geometry dimensions, & viscosity
    rotx(α) = SA{T}[1 0 0; 0 cos(α) -sin(α); 0 sin(α) cos(α)]
    roty(α) = SA{T}[cos(α) 0 sin(α); 0 1 0; -sin(α) 0 cos(α)]
    rotz(α) = SA{T}[cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1]
    function valve(x,t,R,H;β=0)
        a = 3.1f0/13.5f0*R; b = 2.f0/13.5f0*R # scale to new radius R
        xᵢ = rotx(atan(R/2/H))*rotz(β)*x.-SA{T}[0,3R/4,0] # rotate
        # cut for points outside √(x²+y²) = R and above H
        return max(max(abs(xᵢ[1]^2/a^2 + xᵢ[2]^2/b^2- xᵢ[3])-2, √sum(abs2,x[1:2])-R), x[3]-H)
    end
    R = 2L/3.f0; H = 18.0f0/13.5f0*R
    body = AutoBody((x,t)->valve(x,t,R,H;β=0),(x,t)->x.-SA[L,L,L-R/4.f0])
    # body += AutoBody((x,t)->valve(x,t,R,H;β=0),(x,t)->x.-SA{T}[L,L+R/2,L]) # add second leaflet
    # body += AutoBody((x,t)->valve(x,t,R,H;β=-2π/3),(x,t)->x.-SA{T}[L-R,L-R,L]) # add third leaflet
   
    # Return initialized simulation
    Simulation((2L,2L,2L),(0,0,0),L;U,ν=U*L/Re,body,mem,T)
end

# using CUDA
# CUDA.allowscalar(false)
# Define geometry and motion on GPU    
sim = leaflet(2^6)#mem=CUDA.CuArray);
wr = vtkWriter("leaflet"; attrib=custom_attrib)
save!(wr, sim)
close(wr)
# Loop in time
# # record(fig,"jelly.mp4",1:200) do frame
# foreach(1:1) do frame
#     @show frame
#     sim_step!(sim,sim_time(sim)+0.05,verbose=true);
#     write!(wr, sim)
# end
# close(wr)