
using VortexCollisions

include("test_functions.jl")


function run_timings()
    M = 64
    N = 64

    grid = Grid2d(M,N)
    ft   = FourierTransform(grid)
    op   = FokkerPlanckOperator(grid, ft)
    D    = get_gradient(ft)
    Δ⁻¹  = get_inverse_laplacian(ft)

    u = get_field(grid)
    û = get_trans(ft)

    Δ⁻¹û = zeros(û)
    Dû   = [zeros(û), zeros(û)]

    Δ⁻¹u = zeros(u)
    divJ = zeros(u)

    evaluate_function_on_grid(grid, u_test, u)


    frfft!(ft, u, û)
    prfft!(ft, u, û)
    apply_operator!(Δ⁻¹, û, Δ⁻¹û)
    apply_operator!(D, û, Dû)
    irfft!(ft, Δ⁻¹û, Δ⁻¹u)
    fourier_quadrature(û, û, ft.μ, grid)
    collision_operator!(op, u, divJ)


    print(" frfft!:     ")
    @time frfft!(ft, u, û)

    print(" prfft!:     ")
    @time prfft!(ft, u, û)

    print(" Δ⁻¹:        ")
    @time apply_operator!(Δ⁻¹, û, Δ⁻¹û)

    print(" D:          ")
    @time apply_operator!(D, û, Dû)

    print(" irfft!:     ")
    @time irfft!(ft, Δ⁻¹û, Δ⁻¹u)

    print(" quadrature: ")
    @time fourier_quadrature(û, û, ft.μ, grid)

    print(" collisions: ")
    @time collision_operator!(op, u, divJ)

end


run_timings()
