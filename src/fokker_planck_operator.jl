

function mfunc(u::Matrix)
    u
end


function hfunc(u::Matrix, ϕ::Matrix)
    ϕ
end


# g_kij, with k the vector component, and ij the grid index
function wfunc{T}(x::T, x_grid::Matrix{T}, g::Array{T,3})
    local g¹::T
    local g²::T

    Π = zeros(T, size(g,1), size(g,1), size(g,2), size(g,3))

    for j in 1:size(g,3)
        for i in 1:size(g,2)
            g² = 0
            for k in 1:size(g,1)
                g² += g[k,i,j]^2
                Π[k,k,i,j] = 1
            end
            g¹ = sqrt(g²)
            for k in 1:size(Π,1)
                for l in 1:size(Π,2)
                    Π[l,k,i,j] += g[l,i,j] * g[k,i,j] / g²
                    Π[l,k,i,j] /= g¹
                end
            end
        end
    end

    Π
end
