
struct FokkerPlanckOperator{RT,CT,M,N,ℳ,𝒩,HT,MT,WT} <: CollisionOperator{RT,CT,M,N,ℳ,𝒩}
    grid::Grid2d{M,N,RT}
    ft::FourierTransform{ℳ,𝒩,RT,CT}
    D::Vector{Matrix{CT}}
    Δ⁻¹::Matrix{CT}

    hfunc::HT
    mfunc::MT
    wfunc::WT

    function FokkerPlanckOperator{RT,CT,M,N,ℳ,𝒩,HT,MT,WT}(grid::Grid2d{M,N,RT}, ft::FourierTransform{ℳ,𝒩,RT,CT},
                                                           hfunc::HT, mfunc::MT, wfunc::WT) where {RT,CT,M,N,ℳ,𝒩,HT,MT,WT}
        D   = get_gradient(ft)
        Δ⁻¹ = get_inverse_laplacian(ft, sign=-1)
        new(grid, ft, D, Δ⁻¹, hfunc, mfunc, wfunc)
    end
end

function FokkerPlanckOperator{RT,CT,M,N,ℳ,𝒩}(grid::Grid2d{M,N,RT}, ft::FourierTransform{ℳ,𝒩,RT,CT};
                                              hfunc=hfunc_fokker_planck!, mfunc=mfunc_fokker_planck!, wfunc=wfunc_fokker_planck!)
    FokkerPlanckOperator{RT,CT,M,N,ℳ,𝒩, typeof(hfunc), typeof(mfunc), typeof(wfunc)}(grid, ft, hfunc, mfunc, wfunc)
end


@generated function collision_operator!{RT,CT,M,N,ℳ,𝒩}(op::FokkerPlanckOperator{RT,CT,M,N,ℳ,𝒩}, u::Matrix{RT}, divJ::Matrix{RT})

    local ϕ::Matrix{RT}
    local m::Matrix{RT}
    local h::Matrix{RT}

    local Du::Vector{Matrix{RT}}
    local Dh::Vector{Matrix{RT}}
    local J::Vector{Matrix{RT}}
    local g::Vector{Matrix{RT}}
    local w::Matrix{Matrix{RT}}

    local 𝔽::Vector{Matrix{RT}}
    local 𝔻::Matrix{Matrix{RT}}

    local û::Matrix{CT}
    local ϕ̂::Matrix{CT}
    local m̂::Matrix{CT}
    local ĥ::Matrix{CT}

    local Dû::Vector{Matrix{CT}}
    local Dĥ::Vector{Matrix{CT}}
    local Ĵ::Vector{Matrix{CT}}
    local ŵ::Matrix{Matrix{CT}}

    local divĴ::Matrix{CT}


    ϕ  = zeros(RT,M,N)
    m  = zeros(RT,M,N)
    h  = zeros(RT,M,N)

    Du = [zeros(RT,M,N), zeros(RT,M,N)]
    Dh = [zeros(RT,M,N), zeros(RT,M,N)]
    J  = [zeros(RT,M,N), zeros(RT,M,N)]
    g  = [zeros(RT,M,N), zeros(RT,M,N)]

    w = Array{Array{RT,2}}(2,2)

    for l in 1:size(w,2)
        for k in 1:size(w,1)
            w[k,l] = zeros(RT,M,N)
        end
    end

    𝔽 = [zeros(RT,M,N), zeros(RT,M,N)]

    𝔻 = Array{Array{RT,2}}(2,2)

    for l in 1:size(𝔻,2)
        for k in 1:size(𝔻,1)
            𝔻[k,l] = zeros(RT,M,N)
        end
    end

    û  = zeros(CT,ℳ,𝒩)
    ϕ̂  = zeros(CT,ℳ,𝒩)
    m̂  = zeros(CT,ℳ,𝒩)
    ĥ  = zeros(CT,ℳ,𝒩)

    Dû = [zeros(CT,ℳ,𝒩), zeros(CT,ℳ,𝒩)]
    Dĥ = [zeros(CT,ℳ,𝒩), zeros(CT,ℳ,𝒩)]
    Ĵ  = [zeros(CT,ℳ,𝒩), zeros(CT,ℳ,𝒩)]

    ŵ = Array{Array{CT,2}}(2,2)

    for l in 1:size(ŵ,2)
        for k in 1:size(ŵ,1)
            ŵ[k,l] = zeros(CT,ℳ,𝒩)
        end
    end

    divĴ = zeros(CT,ℳ,𝒩)


    quote
        @assert size(u) == size(divJ) == (M,N)

        frfft!(op.ft, u, $û)

        apply_operator!(op.Δ⁻¹, $û, $ϕ̂)

        irfft!(op.ft, $ϕ̂, $ϕ)

        op.mfunc(u, $m)
        op.hfunc(u, $ϕ, $h)

        frfft!(op.ft, $m, $m̂)
        frfft!(op.ft, $h, $ĥ)

        apply_operator!(op.D, $û, $Dû)
        apply_operator!(op.D, $ĥ, $Dĥ)

        for k in 1:length($Dĥ)
            irfft!(op.ft, $Dĥ[k], $Dh[k])
        end


        for j in 1:N
            for i in 1:M

                for k in 1:length($g)
                    $g[k] .= $Dh[k][i,j] .- $Dh[k]
                end

                op.wfunc(op.grid, i, j, $g, $w)

                for l in 1:size($w,2)
                    for k in 1:size($w,1)
                        prfft!(op.ft, $w[k,l], $ŵ[k,l])
                    end
                end

                for k in 1:size($ŵ,1)
                    $𝔽[k][i,j] = 0
                    for l in 1:size($ŵ,2)
                        $𝔽[k][i,j] -= real(fourier_quadrature($ŵ[k,l], $Dû[l], op.ft.μ, op.grid))
                    end
                end

                for l in 1:size($𝔻,2)
                    for k in 1:size($𝔻,1)
                        $𝔻[k,l][i,j] = real(fourier_quadrature($ŵ[k,l], $m̂, op.ft.μ, op.grid))
                    end
                end

            end
        end


        for k in 1:length($Du)
            irfft!(op.ft, $Dû[k], $Du[k])
        end

        for k in 1:length($J)
            for j in 1:size($J[k],2)
                for i in 1:size($J[k],1)
                    $J[k][i,j] = $𝔻[k,1][i,j] * $Du[1][i,j] + $𝔻[k,2][i,j] * $Du[2][i,j] + $𝔽[k][i,j] * $m[i,j]
                end
            end
            prfft!(op.ft, $J[k], $Ĵ[k])
        end

        apply_operator!(op.D, $Ĵ, $divĴ)

        irfft!(op.ft, $divĴ, divJ)

        divJ
    end
end


function mfunc_fokker_planck!(u::Matrix, m::Matrix)
    m .= 1
end


function hfunc_fokker_planck!(u::Matrix, ϕ::Matrix, h::Matrix)
    h .= ϕ
end


@generated function wfunc_fokker_planck!{M,N,T}(grid::Grid2d{M,N,T}, i::Int, j::Int, g::Vector{Matrix{T}}, w::Matrix{Matrix{T}})
    # local g¹::Matrix{T} = zeros(T,M,N)
    local g²::Matrix{T} = zeros(T,M,N)

    quote
        @assert size(w, 1) == length(g)
        @assert size(w, 2) == length(g)
        @assert size(w[1,1], 1) == size(g[1],1) == M
        @assert size(w[1,1], 2) == size(g[1],2) == N

        for l in 1:size(w,2)
            for k in 1:size(w,1)
                w[k,l] .= 0
            end
        end

        $g² .= 0

        for k in 1:size(g,1)
            for j in 1:size($g²,2)
                for i in 1:size($g²,1)
                    $g²[i,j] += g[k][i,j]^2
                end
            end
        end

        # $g¹ .= sqrt.($g²)

        for k in 1:size(w,1)
            for j in 1:size(w[k,k],2)
                for i in 1:size(w[k,k],1)
                    w[k,k][i,j] = $g²[i,j]
                end
            end
        end

        for l in 1:size(w,2)
            for k in 1:size(w,1)
                for j in 1:size(w[k,l],2)
                    for i in 1:size(w[k,l],1)
                        w[k,l][i,j] -= g[k][i,j] * g[l][i,j]
                    end
                end
            end
        end

        w
    end
end
