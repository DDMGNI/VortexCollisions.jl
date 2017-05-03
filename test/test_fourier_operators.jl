
function u_test(x, y)
    exp(cos(x) + cos(2y))
end

function dudx_test(x,y)
    - sin(x) * u_test(x,y)
end

function dudy_test(x,y)
    - 2sin(2y) * u_test(x,y)
end

function testGradient()
    M = 32
    N = 64

    grid = Grid2d(M,N)
    ft   = FourierTransform(Float64, grid)
    D    = get_gradient(ft)

    u  = get_field(grid)
    ux = zeros(u)
    uy = zeros(u)

    for i in 1:M
        for j in 1:N
            u[i,j]  = u_test(grid.x[i], grid.y[j])
            ux[i,j] = dudx_test(grid.x[i], grid.y[j])
            uy[i,j] = dudy_test(grid.x[i], grid.y[j])
        end
    end

    û = get_trans(ft)
    prfft!(ft, u, û)
    Dû = [zeros(û), zeros(û)]
    apply_operator!(D, û, Dû)

    uxfft = zeros(ux)
    uyfft = zeros(uy)

    irfft!(ft, Dû[1], uxfft)
    irfft!(ft, Dû[2], uyfft)

    @test maximum(ux - uxfft) ≈ zero(eltype(u)) atol=1E-12
    @test maximum(uy - uyfft) ≈ zero(eltype(u)) atol=1E-12

end


function f_test(x, y)
    (sin(x)^2 - cos(x) + 4sin(2y)^2 - 4cos(2y)) * u_test(x,y)
end

function testInverseLaplacian()
    M = 32
    N = 64

    grid = Grid2d(M,N)
    ft   = FourierTransform(Float64, grid)
    Δ⁻¹  = get_inverse_laplacian(ft)

    u = get_field(grid)
    f = get_field(grid)

    for i in 1:M
        for j in 1:N
            u[i,j]  = u_test(grid.x[i], grid.y[j])
            f[i,j]  = f_test(grid.x[i], grid.y[j])
        end
    end

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

    @test maximum(Δ⁻¹f_ana - Δ⁻¹f_fft) ≈ zero(eltype(f)) atol=1E-14

end


testGradient()
testInverseLaplacian()
