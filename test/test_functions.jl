
function cos_test(x, y)
    cos(x) * cos(y)
end

function u_test(x, y)
    exp(cos(x) + cos(2y))
end

function dudx_test(x,y)
    - sin(x) * u_test(x,y)
end

function dudy_test(x,y)
    - 2sin(2y) * u_test(x,y)
end


function f_test(x, y)
    (sin(x)^2 - cos(x) + 4sin(2y)^2 - 4cos(2y)) * u_test(x,y)
end


function evaluate_function_on_grid(grid, ufunc, u)
    @assert size(u,1) == length(grid.x)
    @assert size(u,2) == length(grid.y)
    @inbounds for j in 1:size(u,2)
        for i in 1:size(u,1)
            u[i,j]  = ufunc(grid.x[i], grid.y[j])
        end
    end
end
