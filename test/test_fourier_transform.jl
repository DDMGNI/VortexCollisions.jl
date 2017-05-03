
function testFourierTransform()
    M = 4
    N = 5

    grid = Grid2d(M,N)
    ft   = FourierTransform(Float64, grid; mcut=2, ncut=2)

    println(ft)
end


function testRealFourierTransforms()
    M  = 32
    N  = 32

    grid = Grid2d(M,N)
    u    = get_field(grid)
    ufft = get_field(grid)
    ft   = FourierTransform(eltype(u), grid; mcut=4, ncut=4)
    û    = get_trans(ft)

    for j in 1:N
        for i in 1:M
            u[i,j]  = cos(grid.x[i]) * sin(grid.y[j])
        end
    end


    prfft!(ft, u, û)
    irfft!(ft, û, ufft)

    @test maximum(u - ufft) ≈ zero(eltype(u)) atol=1E-15


    frfft!(ft, u, û)
    irfft!(ft, û, ufft)

    @test maximum(u - ufft) ≈ zero(eltype(u)) atol=1E-15

end


testFourierTransform()
testRealFourierTransforms()
