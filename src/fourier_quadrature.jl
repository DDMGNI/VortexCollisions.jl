

function fourier_quadrature{T}(w::Matrix{T}, v::Matrix{T}, u::Matrix{T}, normal)
    @assert size(u) == size(v) == size(w)

    local result::T = 0

    for j in 1:size(u,2)
        for i in 1:size(u,1)
            result += w[i,j] * v[i,j] * u[i,j]
        end
    end

    result * normal
end


function fourier_quadrature{T}(w::Matrix{Matrix{T}}, v::Matrix{T}, u::Matrix{T}, normal)
    local result::Matrix{T} = zeros(T, size(w, 1), size(w, 2))

    for j in 1:size(w,2)
        for i in 1:size(w,1)
            result[i,j] = fourier_quarature(w[i,j], v, u, normal)
        end
    end

    result
end


function fourier_quadrature{T}(w::Matrix{Matrix{T}}, v::Vector{Matrix{T}}, u::Matrix{T}, normal)
    @assert size(w, 2) == length(v)

    local result::Vector{T} = zeros(T, size(w, 1))

    for i in 1:size(w,1)
        for j in 1:size(w,2)
            result[i] += fourier_quarature(w[i,j], v[j], u, normal)
        end
    end

    result
end
