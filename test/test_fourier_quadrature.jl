
function u_test(x, y, σ)
    exp( - ( (x-π)^2 + (y-π)^2 ) / (2σ^2) )
end

function testQuadrature()
    M = 64
    N = 64

    σ = 0.5

    grid = Grid2d(M,N)
    fmp = FourierMultipliers(Float64, grid)

    u = zeros(M,N)

    for i in 1:M
        for j in 1:N
            u[i,j]  = u_test(grid.x[i], grid.y[j], σ)
        end
    end

    û    = rfft(u, (2,1))
    ûint = fourier_quadrature(û, û, fmp.μ, grid.normalisation)

    @test real(ûint) - π*σ^2 ≈ zero(eltype(u)) atol=1E-15
    @test imag(ûint) ≈ zero(eltype(u)) atol=eps()
end


testQuadrature()
