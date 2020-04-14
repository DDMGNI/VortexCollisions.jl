

function trapezoidal_quadrature(w::Matrix{RT}, v::Union{Array{RT,2},SharedArray{RT,2}}, grid::Grid2d{M,N,RT}) where {M,N,RT}
    @assert size(v) == size(w)

    local result::RT = 0

    # @inbounds for j in 1:size(v,2)
    #     for i in 1:size(v,1)
    #         result += w[i,j] * v[i,j]
    #     end
    # end

    @inbounds for i in eachindex(v,w)
          result += w[i] * v[i]
    end

    result * (grid.x2 - grid.x1) / M * (grid.y2 - grid.y1) / N
end


function trapezoidal_quadrature(w::Matrix{Matrix{RT}}, v::Union{Array{RT,2},SharedArray{RT,2}}, grid::Grid2d{M,N,RT}) where {M,N,RT}
    local result::Matrix{CT} = zeros(RT, size(w,1), size(w,2))

    for l in 1:size(w,2)
        for k in 1:size(w,1)
            result[k,l] = trapezoidal_quadrature(w[k,l], v, grid)
        end
    end

    result
end


function trapezoidal_quadrature(w::Matrix{Matrix{RT}}, v::Union{Vector{Array{RT,2}},Vector{SharedArray{RT,2}}}, grid::Grid2d{M,N,RT}) where {M,N,RT}
    @assert size(w,2) == length(v)

    local result::Vector{CT} = zeros(RT, size(w,1))

    for k in 1:size(w,1)
        for l in 1:size(w,2)
            result[k] += trapezoidal_quadrature(w[k,l], v[l], grid)
        end
    end

    result
end
