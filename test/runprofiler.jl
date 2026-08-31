
# call with julia --track-allocation=user test/runprofiler.jl

using VortexCollisions

include("test_functions.jl")

function run_profiler()
    M = 64
    N = 64

    ℳcut = M
    𝒩cut = N

    grid = Grid2d(M, N)
    ft = FourierTransform(grid; ℳcut = ℳcut, 𝒩cut = 𝒩cut)
    op = FokkerPlanckOperator(grid, ft)

    u = get_field(grid)
    divJ = zeros(u)

    evaluate_function_on_grid(grid, u_test, u)

    collision_operator!(op, u, divJ)

    Profile.clear_malloc_data()

    collision_operator!(op, u, divJ)
end

run_profiler()
