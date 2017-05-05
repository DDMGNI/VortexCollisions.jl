
function u_test(x, y, σ)
    exp( - ( (x-π)^2 + (y-π)^2 ) / (2σ^2) )
end

function testTrapezoidalQuadrature()
    M = 64
    N = 64

    σ = 0.5

    grid = Grid2d(M,N)
    u    = get_field(grid)

    for i in 1:M
        for j in 1:N
            u[i,j]  = u_test(grid.x[i], grid.y[j], σ)
        end
    end

    uint = trapezoidal_quadrature(u, u, grid)

    @test uint - π*σ^2 ≈ zero(eltype(u)) atol=1E-14
end


testTrapezoidalQuadrature()
