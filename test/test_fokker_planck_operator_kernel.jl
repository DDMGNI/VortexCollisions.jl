
using VortexCollisions

include("test_functions.jl")

function u_test_sinx4(x, y)
    sin(x)^4
end

function mfunc_one!(u, m)
    m .= 1
end

function hfunc_ϕ!(u, ϕ, h)
    h .= ϕ
end


function testNullspace(utest, mfunc, hfunc)
    M  = 64
    N  = 64

    grid = Grid2d(M,N)
    u    = get_field(grid)
    divJ = get_field(grid)
    ft   = FourierTransform(grid)
    op   = FokkerPlanckOperator(grid, ft; mfunc=mfunc, hfunc=hfunc)

    evaluate_function_on_grid(grid, utest, u)

    collision_operator!(op, u, divJ)

    @test maximum(abs.(divJ)) == 0.0

end


testNullspace(u_test_sinx4, mfunc_one!, hfunc_ϕ!)
