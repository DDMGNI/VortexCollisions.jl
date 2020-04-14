
using AbstractFFTs
using FFTW


struct FourierTransform{ℳ, 𝒩, RT <: Number, CT <: Number}
    ℳcut::Int
    𝒩cut::Int
    m::Vector{RT}
    n::Vector{RT}
    χ::Matrix{RT}
    μ::Matrix{CT}
    forw_plan::FFTW.rFFTWPlan{RT,-1,false,2}
    back_plan::AbstractFFTs.ScaledPlan{CT,FFTW.rFFTWPlan{CT,1,false,2},RT}
end

function FourierTransform(grid::Grid2d{M,N,RT}; ℳcut::Int=0, 𝒩cut::Int=0) where {M,N,RT}
    ℳ = M
    𝒩 = div(N,2)+1

    CT = typeof(complex(one(RT),0))

    m = circshift(collect(-div(ℳ-1,2):div(ℳ,2)), -div(ℳ-1,2))
    n = collect(1:𝒩) .- 1

    χ   = cutoff_frequencies(RT, ℳ, 𝒩, ℳcut, 𝒩cut)
    μ   = multiplicity(CT, ℳ, 𝒩)

    forw_plan = plan_rfft(zeros(RT, M, N), (2,1))
    back_plan = plan_irfft(zeros(CT, ℳ, 𝒩), N, (2,1))

    # println(size(forw_plan), eltype(forw_plan), typeof(forw_plan))
    # println(size(back_plan), eltype(back_plan), typeof(back_plan))

    FourierTransform{ℳ, 𝒩, RT, CT}(ℳcut, 𝒩cut, m, n, χ, μ, forw_plan, back_plan)
end

function Base.show(io::IO, ft::FourierTransform{ℳ,𝒩,RT,CT}) where {ℳ,𝒩,RT,CT}
    print(io, "Fourier Transform for data types (", RT, ",", CT, ") with (", ℳ, ",", 𝒩, ") frequencies.\n")
    print(io, "   Cut off frequencies: ", ft.χ, "\n")
    print(io, "   Multiplicity:        ", ft.μ, "\n")
end


function get_trans(ft::FourierTransform{ℳ,𝒩,RT,CT}) where {ℳ,𝒩,RT,CT}
    zeros(CT, ℳ, 𝒩)
end


"""
Plain (real) FFT.
"""
function prfft!(ft::FourierTransform{ℳ,𝒩,RT,CT}, u::Matrix{RT}, û::Matrix{CT}) where {ℳ,𝒩,RT,CT}
    mul!(û, ft.forw_plan, u)
end


"""
Filtered (real) FFT.
"""
function frfft!(ft::FourierTransform{ℳ,𝒩,RT,CT}, u::Matrix{RT}, û::Matrix{CT}) where {ℳ,𝒩,RT,CT}
    mul!(û, ft.forw_plan, u)
    û .*= ft.χ
end

"""
Inverse (real) FFT.
"""
function irfft!(ft::FourierTransform{ℳ,𝒩,RT,CT}, û::Matrix{CT}, u::Matrix{RT}) where {ℳ,𝒩,RT,CT}
    mul!(u, ft.back_plan, û)
end


function cutoff_frequencies(T, ℳ, 𝒩, ℳcut, 𝒩cut)
    χ = ones(T, ℳ, 𝒩)

    if ℳcut > 0
        for n in 1:𝒩
            for m in ℳcut+1:div(ℳ+1,2)
                χ[m, n] = 0
            end
            for m in ℳcut-1:div(ℳ,2)
                χ[end-m, n] = 0
            end
        end
    end

    if 𝒩cut > 0
        for n in 𝒩cut+1:𝒩
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
