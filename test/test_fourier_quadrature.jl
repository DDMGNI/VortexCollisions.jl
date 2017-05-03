
function u_test(x, y, σ)
    exp( - ( (x-π)^2 + (y-π)^2 ) / (2σ^2) )
end

function testQuadrature()
    M = 64
    N = 64

    σ = 0.5

    grid = Grid2d(M,N)
    ft   = FourierTransform(grid)
    u    = get_field(grid)
    û    = get_trans(ft)

    for i in 1:M
        for j in 1:N
            u[i,j]  = u_test(grid.x[i], grid.y[j], σ)
        end
    end

    prfft!(ft, u, û)
    ûint = fourier_quadrature(û, û, ft.μ, grid)

    @test real(ûint) - π*σ^2 ≈ zero(eltype(u)) atol=1E-15
    @test imag(ûint) ≈ zero(eltype(u)) atol=eps()
end


testQuadrature()
