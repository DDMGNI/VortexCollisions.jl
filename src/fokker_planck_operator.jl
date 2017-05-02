
struct FokkerPlanckOperator{T,M,N}
    grid::Grid2d{M,N}
    g::Array{T,3}
    w::Array{T,4}
end

function FokkerPlanckOperator{M,N}(T, grid::Grid2d{M,N})
    g = zeros(T, 2, M, N)
    w = zeros(T, 2, 2, M, N)

    FokkerPlanckOperator{T,M,N}(grid, g, w)
end


function mfunc(op::FokkerPlanckOperator, u::Matrix)
    u
end


function hfunc(op::FokkerPlanckOperator, u::Matrix, ϕ::Matrix)
    ϕ
end


# g_kij, with k the vector component, and ij the grid index
function wfunc{T}(op::FokkerPlanckOperator{T}, x::T, x_grid::Matrix{T}, g::Array{T,3}, w::Array{T,4})
    local g¹::T
    local g²::T

    @assert size(w, 1) == size(g,1)
    @assert size(w, 2) == size(g,1)
    @assert size(w, 3) == size(g,2)
    @assert size(w, 4) == size(g,3)

    for j in 1:size(g,3)
        for i in 1:size(g,2)
            g² = 0
            for k in 1:size(g,1)
                g² += g[k,i,j]^2
                w[k,k,i,j] = 1
            end
            g¹ = sqrt(g²)
            for k in 1:size(w,1)
                for l in 1:size(w,2)
                    w[l,k,i,j] += g[l,i,j] * g[k,i,j] / g²
                    w[l,k,i,j] /= g¹
                end
            end
        end
    end

    w
end
