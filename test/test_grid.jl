
function testGrid()
    grid = Grid2d(7,7)
    @assert frequencies(1, 1, grid) == ( 0,0)
    @assert frequencies(2, 3, grid) == ( 1,2)
    @assert frequencies(5, 4, grid) == (-3,3)
end


testGrid()
