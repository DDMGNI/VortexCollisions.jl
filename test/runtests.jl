
using VortexCollisions


function testFourierMultipliers()
    grid = Grid2d(7,6)
    frequencies(1, 1, grid)
    frequencies(2, 3, grid)
    fmp = FourierMultipliers(Float64, 6, 6, 2, 2)
    println(fmp)
end


function testCollisionMultipliers()
    grid = Grid2d(7,6)
    op = FokkerPlanckOperator(Float64, grid)
end


testFourierMultipliers()
testCollisionMultipliers()
