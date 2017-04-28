
immutable Grid{M,N}

end


immutable FourierMultipliers{RT, CT}
    χ::Matrix{RT}
    μ::Matrix{RT}
    Δ⁻¹::Matrix{RT}
    D::Vector{Matrix{CT}}
    normal::RT
end

function FourierMultipliers(M, N, mcut, ncut)
    ℳ = M
    𝒩 = div(N,2)+1

    RT = Float64
    CT = Complex128

    χ   = cutoff_frequencies(RT, ℳ, 𝒩, mcut, ncut)
    μ   = multiplicity(RT, ℳ, 𝒩)
    Δ⁻¹ = laplacian(RT, M, N, ℳ, 𝒩)
    D   = gradient(CT, M, N, ℳ, 𝒩)

    normal = 4π^2 / M^2 / N^2

    # println(χ)
    # println(μ)
    # println(Δ⁻¹)
    # println(D)

    FourierMultipliers{RT,CT}(χ, μ, Δ⁻¹, D, normal)
end


function cutoff_frequencies(T, ℳ, 𝒩, mcut, ncut)
    χ = ones(T, ℳ, 𝒩)

    for n in 1:𝒩
        for m in mcut+1:div(ℳ+1,2)
            χ[m, n] = 0
        end
        for m in mcut-1:div(ℳ,2)
            χ[end-m, n] = 0
        end
    end

    for n in ncut+1:𝒩
        for m in 1:ℳ
            χ[m,n] = 0
        end
    end

    χ
end


function multiplicity(T, ℳ, 𝒩)
    μ = 2ones(T, ℳ, 𝒩)
    μ[:,1] .= 1
    μ
end


function laplacian(T, M, N, ℳ, 𝒩)
    m = circshift(collect(-div(M-1,2):div(M,2)), -div(M-1,2))
    n = collect(0:div(N,2))

    Δ⁻¹ = zeros(T, ℳ, 𝒩)

    for s in 2:𝒩
        for r in 2:ℳ
            Δ⁻¹[r,s] = - 1 / (m[r]^2 + n[s]^2)
        end
    end

    Δ⁻¹
end


function gradient(T, M, N, ℳ, 𝒩)
    m = circshift(collect(-div(M-1,2):div(M,2)), -div(M-1,2))
    n = collect(0:div(N,2))

    D₁ = zeros(T, ℳ, 𝒩)
    D₂ = zeros(T, ℳ, 𝒩)

    for s in 2:𝒩
        for r in 2:ℳ
            D₁[r,s] = 1im * m[r]
            D₂[r,s] = 1im * n[s]
        end
    end

    [D₁, D₂]
end


@generated function frequencies{M, N, T <: Integer}(r::T, s::T, grid::Grid{M,N})
    m = circshift(collect(-div(M-1,2):div(M,2)), -div(M-1,2))
    n = collect(0:div(N,2))

    # println(m)
    # println(n)

    quote
        return $m[r], $n[s]
    end
end
