


"""
Apply fourier operator v ← op * u
"""
function apply_operator!{CT}(op::Matrix{CT}, u::Matrix{CT}, v::Matrix{CT})
    v .= op .* u
end

"""
Apply fourier operator v[i] ← op[i] * u for i in 1:length(op)
"""
function apply_operator!{CT}(op::Vector{Matrix{CT}}, u::Matrix{CT}, v::Vector{Matrix{CT}})
    for k in 1:length(op)
        v[k] .= op[k] .* u
    end
end


function get_gradient{ℳ,𝒩,RT,CT}(ft::FourierTransform{ℳ,𝒩,RT,CT})
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


function get_inverse_laplacian{ℳ,𝒩,RT,CT}(ft::FourierTransform{ℳ,𝒩,RT,CT})
    Δ⁻¹ = zeros(CT, ℳ, 𝒩)

    for s in 1:𝒩
        for r in 1:ℳ
            if !(r == 1 && s == 1)
                Δ⁻¹[r,s] = - 1 / (ft.m[r]^2 + ft.n[s]^2)
            end
        end
    end

    Δ⁻¹
end
