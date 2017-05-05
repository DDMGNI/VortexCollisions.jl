

function fourier_quadrature{M,N,RT,CT}(w::Matrix{CT}, v::Union{Array{CT,2},SharedArray{CT,2}}, u::Union{Array{CT,2},SharedArray{CT,2}}, grid::Grid2d{M,N,RT})
    @assert size(u) == size(v) == size(w)

    local result::CT = 0

    @inbounds for j in 1:size(u,2)
        for i in 1:size(u,1)
            result += w[i,j] * conj(v[i,j]) * u[i,j]
        end
    end

    result * grid.normalisation
end


function fourier_quadrature{M,N,RT,CT}(w::Matrix{Matrix{CT}}, v::Union{Array{CT,2},SharedArray{CT,2}}, u::Union{Array{CT,2},SharedArray{CT,2}}, grid::Grid2d{M,N,RT})
    local result::Matrix{CT} = zeros(CT, size(w,1), size(w,2))

    for l in 1:size(w,2)
        for k in 1:size(w,1)
            result[k,l] = fourier_quarature(w[k,l], v, u, grid)
        end
    end

    result
end


function fourier_quadrature{M,N,RT,CT}(w::Matrix{Matrix{CT}}, v::Union{Vector{Array{CT,2}},Vector{SharedArray{CT,2}}}, u::Union{Array{CT,2},SharedArray{CT,2}}, grid::Grid2d{M,N,RT})
    @assert size(w,2) == length(v)

    local result::Vector{CT} = zeros(CT, size(w,1))

    for k in 1:size(w,1)
        for l in 1:size(w,2)
            result[k] += fourier_quarature(w[k,l], v[l], u, grid)
        end
    end

    result
end


# function fourier_quadrature{RT,CT}(w::Matrix{CT}, v::Union{Array{CT,3},SharedArray{CT,3}}, l, u::Union{Array{CT,2},SharedArray{CT,2}}, grid::Grid2d{RT})
#     @assert size(u) == (size(v,1), size(v,2)) == size(w)
#     @assert l ≤ size(v,3)
#
#     local result::CT = 0
#
#     @inbounds for j in 1:size(u,2)
#         for i in 1:size(u,1)
#             result += w[i,j] * conj(v[i,j,l]) * u[i,j]
#         end
#     end
#
#     result * grid.normalisation
# end
