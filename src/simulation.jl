
using HDF5

"""
Initialise u₀ from uinit(x,y) function and run simulation.
"""
function run_simulation(op::CollisionOperator, nt::Int, uinit::Function, output::HDF5File)


end


"""
Initialise u₀ from HDF5 file and run simulation.
"""
function run_simulation(op::CollisionOperator, nt::Int, input::HDF5File, output::HDF5File)


end


"""
Run simulation for nt time steps starting from u₀ and write output to HDF5 file.
"""
function run_simulation{RT,CT,M,N,ℳ,𝒩}(op::CollisionOperator{RT,CT,M,N,ℳ,𝒩}, nt::Int, u₀::Matrix{RT}, output::HDF5File)



end
