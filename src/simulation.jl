
using HDF5


"""
Initialise u₀ from uinit(x,y) function and run simulation.
"""
function run_simulation(op::CollisionOperator, nt::Int, Δt::Number, uinit::Function, output::String)
    u₀ = get_field(op.grid)

    for j in 1:size(u₀,2)
        for i in 1:size(u₀,1)
             u₀[i,j]  = uinit(op.grid.x[i], op.grid.y[j])
        end
    end

    run_simulation(op, nt, Δt, u₀, output)
end


"""
Initialise u₀ from HDF5 file <input> and run simulation.
"""
function run_simulation(op::CollisionOperator, nt::Int, Δt::Number, input::String, output::String)
    u₀ = get_field(op.grid)

    h5open(input, "r") do h5
        @assert size(h5.ω,1) == size(op.grid,1)
        @assert size(h5.ω,2) == size(op.grid,2)
        u₀ .= h5.ω[:,:,end]
    end

    run_simulation(op, nt, Δt, u₀, output)
end


"""
Run simulation for nt time steps Δt starting from u₀ and write output to HDF5 file <output>.
"""
function run_simulation{RT,CT,M,N}(op::CollisionOperator{RT,CT,M,N}, nt::Int, Δt::RT, u₀::Matrix{RT}, output::String)

    # create HDF5 output file
    h5 = h5open(output, "w")

    h5ϕ = d_create(h5, "ϕ", datatype(RT), dataspace(M,N,nt+1), "chunk", (M,N,1))
    h5ω = d_create(h5, "ω", datatype(RT), dataspace(M,N,nt+1), "chunk", (M,N,1))


    # write initial conditions to HDF5 file
    write_solution_to_hdf5(op, u₀, 1, h5ϕ, h5ω)


    # run for nt time steps
    u₁ = zeros(u₀)

    for n in 1:nt
        timestep!(op, u₀, u₁, Δt)
        write_solution_to_hdf5(op, u₁, n+1, h5ϕ, h5ω)
        u₀ .= u₁
    end


    # close HDF5 output file
    close(h5)

end


@generated function write_solution_to_hdf5{RT,CT,M,N,ℳ,𝒩}(op::CollisionOperator{RT,CT,M,N,ℳ,𝒩}, u::Matrix{RT}, n::Int, h5ϕ::HDF5Dataset, h5ω::HDF5Dataset)
    local û::Matrix{CT} = zeros(CT,ℳ,𝒩)
    local ϕ̂::Matrix{CT} = zeros(CT,ℳ,𝒩)
    local ϕ::Matrix{RT} = zeros(RT,M,N)

    quote
        prfft!(op.ft, u, $û)
        $ϕ̂ .= op.Δ⁻¹ .* $û
        irfft!(op.ft, $ϕ̂, $ϕ)

        h5ϕ[:,:,n] = $ϕ
        h5ω[:,:,n] = u
    end
end
