
using VortexCollisions


function u_test(x, y)

end


function run_timings()
    M = 32
    N = 64

    grid = Grid2d(M,N)
    ft   = FourierTransform(Float64, grid)
    D    = get_gradient(ft)
    Δ⁻¹  = get_inverse_laplacian(ft)

    u = get_field(grid)
    û = get_trans(ft)

    Δ⁻¹û = zeros(û)
    Dû   = [zeros(û), zeros(û)]

    Δ⁻¹u = zeros(u)

    for i in 1:M
        for j in 1:N
            u[i,j]  = exp(cos(grid.x[i]) + cos(grid.y[j]))
        end
    end


    frfft!(ft, u, û)
    prfft!(ft, u, û)
    apply_operator!(Δ⁻¹, û, Δ⁻¹û)
    apply_operator!(D, û, Dû)
    irfft!(ft, Δ⁻¹û, Δ⁻¹u)
    fourier_quadrature(û, û, ft.μ, grid)


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


end


run_timings()
