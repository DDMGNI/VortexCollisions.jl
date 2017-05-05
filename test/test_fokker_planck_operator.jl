
function u_test1(x, y)
    cos(x) * cos(y)
end

function u_test2(x, y)
    cos(x) - cos(y)
end

function mfunc_one!(u, m)
    m .= 1
end

function hfunc_u!(u, ϕ, h)
    h .= u
end

function hfunc_ϕ!(u, ϕ, h)
    h .= ϕ
end


# function testFokkerPlanckOperator()
#     M  = 64
#     N  = 64
#
#     grid = Grid2d(M,N)
#     u    = get_field(grid)
#     divJ = get_field(grid)
#     ft   = FourierTransform(grid)
#     op   = FokkerPlanckOperator(grid, ft)
#
#     evaluate_function_on_grid(grid, u_test, u)
#
#     collision_operator!(op, u, divJ)
#
# end


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

    # println(maximum(abs.(divJ)))
    @test maximum(abs.(divJ)) ≈ zero(eltype(divJ)) atol=1E-9
    # @test maximum(abs.(divJ)) ≈ zero(eltype(divJ)) atol=1E-10

end


# testFokkerPlanckOperator()
testNullspace(u_test1, mfunc_one!, hfunc_ϕ!)
testNullspace(u_test2, mfunc_one!, hfunc_ϕ!)
testNullspace(u_test,  mfunc_one!, hfunc_u!)
