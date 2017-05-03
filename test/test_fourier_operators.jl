
function testGradient()
    M = 32
    N = 64

    grid = Grid2d(M,N)
    ft   = FourierTransform(grid)
    D    = get_gradient(ft)

    u  = get_field(grid)
    ux = zeros(u)
    uy = zeros(u)

    evaluate_function_on_grid(grid, u_test, u)
    evaluate_function_on_grid(grid, dudx_test, ux)
    evaluate_function_on_grid(grid, dudy_test, uy)

    û = get_trans(ft)
    prfft!(ft, u, û)
    Dû = [zeros(û), zeros(û)]
    apply_operator!(D, û, Dû)

    uxfft = zeros(ux)
    uyfft = zeros(uy)

    irfft!(ft, Dû[1], uxfft)
    irfft!(ft, Dû[2], uyfft)

    @test maximum(abs.(ux - uxfft)) ≈ zero(eltype(u)) atol=1E-12
    @test maximum(abs.(uy - uyfft)) ≈ zero(eltype(u)) atol=1E-12

end


function testInverseLaplacian()
    M = 32
    N = 64

    grid = Grid2d(M,N)
    ft   = FourierTransform(grid)
    Δ⁻¹  = get_inverse_laplacian(ft)

    u = get_field(grid)
    f = get_field(grid)

    evaluate_function_on_grid(grid, u_test, u)
    evaluate_function_on_grid(grid, f_test, f)

    uhat = get_trans(ft)
    fhat = get_trans(ft)

    prfft!(ft, u, uhat)
    prfft!(ft, f, fhat)

    Δ⁻¹fhat = zeros(fhat)

    apply_operator!(Δ⁻¹, fhat, Δ⁻¹fhat)

    u₀hat = zeros(uhat)
    u₀hat[1,1] = uhat[1,1]

    Δ⁻¹f_fft = zeros(u)
    u₀fft    = zeros(u)

    irfft!(ft, Δ⁻¹fhat, Δ⁻¹f_fft)
    irfft!(ft, u₀hat, u₀fft)

    Δ⁻¹f_ana = u .- u₀fft

    @test maximum(abs.(Δ⁻¹f_ana - Δ⁻¹f_fft)) ≈ zero(eltype(f)) atol=1E-14

end


testGradient()
testInverseLaplacian()
