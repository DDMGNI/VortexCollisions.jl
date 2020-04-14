
"""
Apply fourier operator v ← op * u
"""
function apply_operator!(op::Matrix{CT}, u::Union{Matrix{CT},SharedArray{CT,2}}, v::Union{Matrix{CT},SharedArray{CT,2}}) where {CT}
    v .= op .* u
end

"""
Apply fourier operator v[i] ← op[i] * u for i in 1:length(op)
"""
function apply_operator!(op::Vector{Matrix{CT}}, u::Union{Matrix{CT},SharedArray{CT,2}}, v::Union{Vector{Matrix{CT}},Vector{SharedArray{CT,2}}}) where {CT}
    for k in 1:length(op)
        v[k] .= op[k] .* u
    end
end

"""
Apply fourier operator v ← sum_i op[i] * u[i]
"""
function apply_operator!(op::Vector{Matrix{CT}}, u::Union{Vector{Matrix{CT}},Vector{SharedArray{CT,2}}}, v::Union{Matrix{CT},SharedArray{CT,2}}) where {CT}
    @assert length(op) == length(u)
    fill!(v, 0)
    for k in 1:length(op)
        v .+= op[k] .* u[k]
    end
end

"""
Returns gradient operator
"""
function get_gradient(ft::FourierTransform{ℳ,𝒩,RT,CT}) where {ℳ,𝒩,RT,CT}
    D₁ = zeros(CT, ℳ, 𝒩)
    D₂ = zeros(CT, ℳ, 𝒩)

    for s in 1:𝒩
        for r in 1:ℳ
            D₁[r,s] = 1im * ft.m[r]
            D₂[r,s] = 1im * ft.n[s]
        end
    end

    [D₁, D₂]
end

"""
Returns inverse Laplacian operator
"""
function get_inverse_laplacian(ft::FourierTransform{ℳ,𝒩,RT,CT}; sign=+1) where {ℳ,𝒩,RT,CT}
    Δ⁻¹ = zeros(CT, ℳ, 𝒩)

    for s in 1:𝒩
        for r in 1:ℳ
            if !(r == 1 && s == 1)
                Δ⁻¹[r,s] = - sign / (ft.m[r]^2 + ft.n[s]^2)
            end
        end
    end

    Δ⁻¹
end
