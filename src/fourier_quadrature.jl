

function fourier_quadrature{RT,CT}(w::Matrix{CT}, v::Matrix{CT}, u::Matrix{CT}, grid::Grid2d{RT})
    @assert size(u) == size(v) == size(w)

    local result::CT = 0

    for j in 1:size(u,2)
        for i in 1:size(u,1)
            result += w[i,j] * v[i,j] * u[i,j]
        end
    end

    result * grid.normalisation
end


function fourier_quadrature{RT,CT}(w::Matrix{Matrix{CT}}, v::Matrix{CT}, u::Matrix{CT}, grid::Grid2d{RT})
    local result::Matrix{CT} = zeros(CT, size(w,1), size(w,2))

    for j in 1:size(w,2)
        for i in 1:size(w,1)
            result[i,j] = fourier_quarature(w[i,j], v, u, grid)
        end
    end

    result
end


function fourier_quadrature{RT,CT}(w::Matrix{Matrix{CT}}, v::Vector{Matrix{CT}}, u::Matrix{CT}, grid::Grid2d{RT})
    @assert size(w,2) == length(v)

    local result::Vector{CT} = zeros(CT, size(w,1))

    for i in 1:size(w,1)
        for j in 1:size(w,2)
            result[i] += fourier_quarature(w[i,j], v[j], u, grid)
        end
    end

    result
end
