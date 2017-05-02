
immutable FourierMultipliers{RT,CT}
    ℳ
    𝒩
    χ::Matrix{RT}
    μ::Matrix{RT}
    Δ⁻¹::Matrix{RT}
    D::Vector{Matrix{CT}}
    normal::RT
end

function FourierMultipliers{M,N}(RT, grid::Grid2d{M,N}; mcut=NaN, ncut=NaN)
    ℳ = M
    𝒩 = div(N,2)+1

    CT = typeof(complex(one(RT),0))

    χ   = cutoff_frequencies(RT, ℳ, 𝒩, mcut, ncut)
    μ   = multiplicity(RT, ℳ, 𝒩)
    Δ⁻¹ = inverse_laplacian(RT, ℳ, 𝒩)
    D   = gradient(CT, ℳ, 𝒩)

    normal = 4π^2 / M^2 / N^2

    FourierMultipliers{RT,CT}(ℳ, 𝒩, χ, μ, Δ⁻¹, D, normal)
end

function Base.show{RT,CT}(io::IO, fmp::FourierMultipliers{RT,CT})
    print(io, "Fourier Multipliers for data types (", RT, ",", CT, ") with (", fmp.ℳ, ",", fmp.𝒩, ") frequencies.\n")
    print(io, "   Cut off frequencies: ", fmp.χ, "\n")
    print(io, "   Multiplicity:        ", fmp.μ, "\n")
    print(io, "   Inverse Laplacian:   ", fmp.Δ⁻¹, "\n")
    print(io, "   Gradient:            ", fmp.D, "\n")
end


function cutoff_frequencies(T, ℳ, 𝒩, mcut, ncut)
    χ = ones(T, ℳ, 𝒩)

    if !isnan(mcut)
        for n in 1:𝒩
            for m in mcut+1:div(ℳ+1,2)
                χ[m, n] = 0
            end
            for m in mcut-1:div(ℳ,2)
                χ[end-m, n] = 0
            end
        end
    end

    if !isnan(ncut)
        for n in ncut+1:𝒩
            for m in 1:ℳ
                χ[m,n] = 0
            end
        end
    end

    χ
end


function multiplicity(T, ℳ, 𝒩)
    μ = 2ones(T, ℳ, 𝒩)
    μ[:,1] .= 1
    μ
end


function inverse_laplacian(T, ℳ, 𝒩)
    m = circshift(collect(-div(ℳ-1,2):div(ℳ,2)), -div(ℳ-1,2))
    n = collect(1:𝒩)-1

    Δ⁻¹ = zeros(T, ℳ, 𝒩)

    for s in 1:𝒩
        for r in 1:ℳ
            if !(r == 1 && s == 1)
                Δ⁻¹[r,s] = - 1 / (m[r]^2 + n[s]^2)
            end
        end
    end

    Δ⁻¹
end


function gradient(T, ℳ, 𝒩)
    m = circshift(collect(-div(ℳ-1,2):div(ℳ,2)), -div(ℳ-1,2))
    n = collect(1:𝒩)-1

    D₁ = zeros(T, ℳ, 𝒩)
    D₂ = zeros(T, ℳ, 𝒩)

    for s in 1:𝒩
        for r in 1:ℳ
            D₁[r,s] = 1im * m[r]
            D₂[r,s] = 1im * n[s]
        end
    end

    [D₁, D₂]
end
