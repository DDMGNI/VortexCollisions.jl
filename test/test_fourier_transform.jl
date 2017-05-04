
function testFourierTransform()
    M = 4
    N = 5

    grid = Grid2d(M,N)
    ft   = FourierTransform(grid; ℳcut=2, 𝒩cut=2)

    println(ft)
end


function testRealFourierTransforms()
    M  = 32
    N  = 32

    grid = Grid2d(M,N)
    u    = get_field(grid)
    ufft = get_field(grid)
    ft   = FourierTransform(grid; ℳcut=4, 𝒩cut=4)
    û    = get_trans(ft)

    evaluate_function_on_grid(grid, cos_test, u)


    prfft!(ft, u, û)
    irfft!(ft, û, ufft)

    @test maximum(u - ufft) ≈ zero(eltype(u)) atol=1E-15


    frfft!(ft, u, û)
    irfft!(ft, û, ufft)

    @test maximum(u - ufft) ≈ zero(eltype(u)) atol=1E-15

end


testFourierTransform()
testRealFourierTransforms()
