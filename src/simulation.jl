
using HDF5
using ProgressMeter

@generated function write_solution_to_hdf5(
        op::CollisionOperator{M, N, ℳ, 𝒩, RT, CT}, u::Matrix{RT}, n::Int, h5ϕ::HDF5.Dataset,
        h5ω::HDF5.Dataset, h5h::HDF5.Dataset, h5m::HDF5.Dataset) where {M, N, ℳ, 𝒩, RT, CT}
    local û::Matrix{CT} = zeros(CT, ℳ, 𝒩)
    local ϕ̂::Matrix{CT} = zeros(CT, ℳ, 𝒩)
    local ϕ::Matrix{RT} = zeros(RT, M, N)

    quote
        prfft!(op.ft, u, $û)
        $ϕ̂ .= op.Δ⁻¹ .* $û
        irfft!(op.ft, $ϕ̂, $ϕ)

        h5ϕ[:, :, n] = $ϕ
        h5ω[:, :, n] = u
        h5h[:, :, n] = op.h
        h5m[:, :, n] = op.m
    end
end

"""
Run simulation for nt time steps Δt starting from u₀ and write output to HDF5 file <output>.
"""
function run_simulation(op::CollisionOperator{M, N, ℳ, 𝒩, RT, CT}, nt::Int, Δt::RT,
        u₀::Matrix{RT}, output::String, nsave::Int = 1) where {M, N, ℳ, 𝒩, RT, CT}

    # create HDF5 output file
    h5 = h5open(output, "w")

    h5ϕ = create_dataset(h5, "ϕ", RT, ((M, N, nt+1), (M, N, -1)), chunk = (M, N, 1))
    h5ω = create_dataset(h5, "ω", RT, ((M, N, nt+1), (M, N, -1)), chunk = (M, N, 1))

    h5h = create_dataset(h5, "h", RT, ((M, N, nt+1), (M, N, -1)), chunk = (M, N, 1))
    h5m = create_dataset(h5, "m", RT, ((M, N, nt+1), (M, N, -1)), chunk = (M, N, 1))

    # save grid and time steps to HDF5
    h5["x"] = op.grid.x
    h5["y"] = op.grid.y

    # call collision operator in order to initialize h and m
    collision_operator!(op, u₀, zero(u₀))

    # write initial conditions to HDF5 file
    write_solution_to_hdf5(op, u₀, 1, h5ϕ, h5ω, h5h, h5m)

    # run for nt time steps
    u₁ = zero(u₀)

    @showprogress 1 for n in 1:nt
        timestep!(op, u₀, u₁, Δt)
        mod(n, nsave) == 0 || n == nt ?
        write_solution_to_hdf5(op, u₁, n+1, h5ϕ, h5ω, h5h, h5m) : nothing
        u₀ .= u₁

        mod(n, 10*nsave) == 0 ? flush(h5) : nothing
    end

    # close HDF5 output file
    close(h5)
end

"""
Initialise u₀ from uinit(x,y) function and run simulation.
"""
function run_simulation(op::CollisionOperator, nt::Int, Δt::Number,
        uinit::Function, output::String, nsave::Int = 1)
    u₀ = get_field(op.grid)

    for j in 1:size(u₀, 2)
        for i in 1:size(u₀, 1)
            u₀[i, j] = uinit(op.grid.x[i], op.grid.y[j])
        end
    end

    run_simulation(op, nt, Δt, u₀, output, nsave)
end

"""
Initialise u₀ from HDF5 file <input> and run simulation.
"""
function run_simulation(op::CollisionOperator, nt::Int, Δt::Number,
        input::String, output::String, nsave::Int = 1)
    u₀ = get_field(op.grid)

    h5open(input, "r") do h5
        @assert size(h5.ω, 1) == size(op.grid, 1)
        @assert size(h5.ω, 2) == size(op.grid, 2)
        u₀ .= view(h5["ω"][:, :, end], :, :, 1)
    end

    run_simulation(op, nt, Δt, u₀, output, nsave)
end
