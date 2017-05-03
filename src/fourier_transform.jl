
struct FourierTransform{ℳ, 𝒩, RT <: Number, CT <: Number}
    m::Vector{RT}
    n::Vector{RT}

    χ::Matrix{RT}
    μ::Matrix{CT}
end

function FourierTransform{M,N}(RT, grid::Grid2d{M,N}; mcut=NaN, ncut=NaN)
    ℳ = M
    𝒩 = div(N,2)+1

    CT = typeof(complex(one(RT),0))

    m = circshift(collect(-div(ℳ-1,2):div(ℳ,2)), -div(ℳ-1,2))
    n = collect(1:𝒩)-1

    χ   = cutoff_frequencies(RT, ℳ, 𝒩, mcut, ncut)
    μ   = multiplicity(CT, ℳ, 𝒩)

    FourierTransform{ℳ, 𝒩, RT, CT}(m, n, χ, μ)
end

function Base.show{ℳ,𝒩,RT,CT}(io::IO, ft::FourierTransform{ℳ,𝒩,RT,CT})
    print(io, "Fourier Transform for data types (", RT, ",", CT, ") with (", ℳ, ",", 𝒩, ") frequencies.\n")
    print(io, "   Cut off frequencies: ", ft.χ, "\n")
    print(io, "   Multiplicity:        ", ft.μ, "\n")
end


function get_trans{ℳ,𝒩,RT,CT}(ft::FourierTransform{ℳ,𝒩,RT,CT})
    zeros(CT, ℳ, 𝒩)
end


"""
Plain (real) FFT.
"""
function prfft!{ℳ,𝒩,RT,CT}(ft::FourierTransform{ℳ,𝒩,RT,CT}, u::Matrix{RT}, û::Matrix{CT})
    û .= rfft(u, (2,1))
end


"""
Filtered (real) FFT.
"""
function frfft!{ℳ,𝒩,RT,CT}(ft::FourierTransform{ℳ,𝒩,RT,CT}, u::Matrix{RT}, û::Matrix{CT})
    û .= ft.χ .* rfft(u, (2,1))
end

"""
Inverse (real) FFT.
"""
function irfft!{ℳ,𝒩,RT,CT}(ft::FourierTransform{ℳ,𝒩,RT,CT}, û::Matrix{CT}, u::Matrix{RT})
    u .= irfft(û, size(u,2), (2,1))
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
