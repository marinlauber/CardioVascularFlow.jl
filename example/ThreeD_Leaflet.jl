using WaterLily
using StaticArrays
using WriteVTK
include("../src/utils.jl")

function leaflet(L=2^5;Re=5e2,mem=Array,U=1,AR=2)
    # Define simulation size, geometry dimensions, & viscosity
    rotx(α) = SA[1 0 0; 0 cos(α) -sin(α); 0 sin(α) cos(α)]
    roty(α) = SA[cos(α) 0 sin(α); 0 1 0; -sin(α) 0 cos(α)]
    rotz(α) = SA[cos(α) -sin(α) 0; sin(α) cos(α) 0; 0 0 1]
    function valve(x,t,R,H)
        a = 3.1/13.5*R; b = 2.0/13.5*R
        xᵢ = rotx(atan(R/2H))*x
        # cut for points outside √(x²+y²) = R
        return max(abs(xᵢ[1]^2/a^2+ xᵢ[2]^2/b^2- xᵢ[3])-2, √sum(abs2,x[1:2])-R)
    end
    R = 2L/3; H = 19.15/13.5*R
    body = AutoBody((x,t)->valve(x,t,R,H),(x,t)->x.-L)
   
    # Return initialized simulation
    Simulation((2L,2L,2L),(0,0,0),L;U,ν=U*L/Re,body,mem,T=Float32)
end

# using CUDA
# CUDA.allowscalar(false)
# Define geometry and motion on GPU    
sim = leaflet(2^6)#mem=CUDA.CuArray);
wr = vtkWriter("leaflet"; attrib=custom_attrib)
write!(wr, sim)
close(wr)
# Loop in time
# # record(fig,"jelly.mp4",1:200) do frame
# foreach(1:1) do frame
#     @show frame
#     sim_step!(sim,sim_time(sim)+0.05,verbose=true);
#     write!(wr, sim)
# end
# close(wr)