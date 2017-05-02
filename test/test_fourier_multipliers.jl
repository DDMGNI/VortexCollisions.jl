
function testFourierMultipliers()
    grid = Grid2d(7,6)
    frequencies(1, 1, grid)
    frequencies(2, 3, grid)
    fmp = FourierMultipliers(Float64, grid, mcut=2, ncut=2)
    println(fmp)
end


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
    fmp = FourierMultipliers(Float64, grid)

    u  = zeros(M,N)
    ux = zeros(M,N)
    uy = zeros(M,N)

    for i in 1:M
        for j in 1:N
            u[i,j]  = u_test(grid.x[i], grid.y[j])
            ux[i,j] = dudx_test(grid.x[i], grid.y[j])
            uy[i,j] = dudy_test(grid.x[i], grid.y[j])
        end
    end

    uhat  = rfft(u, (2,1))
    uxhat = fmp.D[1] .* uhat
    uyhat = fmp.D[2] .* uhat

    uxfft = irfft(uxhat, N, (2,1))
    uyfft = irfft(uyhat, N, (2,1))

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
    fmp = FourierMultipliers(Float64, grid)

    u  = zeros(M,N)
    f  = zeros(M,N)

    for i in 1:M
        for j in 1:N
            u[i,j]  = u_test(grid.x[i], grid.y[j])
            f[i,j]  = f_test(grid.x[i], grid.y[j])
        end
    end

    uhat  = rfft(u, (2,1))
    fhat  = rfft(f, (2,1))

    Δ⁻¹fhat = fmp.Δ⁻¹ .* fhat

    u₀hat = zeros(uhat)
    u₀hat[1,1] = uhat[1,1]

    Δ⁻¹f_fft = irfft(Δ⁻¹fhat, N, (2,1))
    Δ⁻¹f_ana = u .- irfft(u₀hat, N, (2,1))

    @test maximum(Δ⁻¹f_ana - Δ⁻¹f_fft) ≈ zero(eltype(f)) atol=1E-14

end


testFourierMultipliers()
testGradient()
testInverseLaplacian()
