
using VortexCollisions


function test()
    grid = Grid{7,6}()
    frequencies(1,1,grid)
    frequencies(2,3,grid)
    FourierMultipliers(6, 6, 2, 2)
end


test()
