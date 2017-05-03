
function testFokkerPlanckOperator()
    M  = 64
    N  = 64

    grid = Grid2d(M,N)
    u    = get_field(grid)
    divJ = get_field(grid)
    ft   = FourierTransform(grid; mcut=31, ncut=31)
    op   = FokkerPlanckOperator(grid, ft)

    evaluate_function_on_grid(grid, u_test, u)

    collision_operator!(op, u, divJ)

end


testFokkerPlanckOperator()
