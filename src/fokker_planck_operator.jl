
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
    # local ĝ::Vector{Matrix{CT}}
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
    # ĝ  = [zeros(CT,ℳ,𝒩), zeros(CT,ℳ,𝒩)]

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
                    # prfft!(op.ft, $g[k], $ĝ[k])
                end

                op.wfunc($g, $w)

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

                for l in 1:size($ŵ,2)
                    for k in 1:size($ŵ,1)
                        $𝔻[k,l][i,j] = real(fourier_quadrature($ŵ[k,l], $m̂, op.ft.μ, op.grid))
                    end
                end

            end
        end


        for k in 1:length($Du)
            irfft!(op.ft, $Dû[k], $Du[k])
        end

        for k in 1:length($J)
            for j in 1:N
                for i in 1:M
                    $J[k][i,j] = $𝔻[k,1][i,j] * $Du[1][i,j] + $𝔻[k,2][i,j] * $Du[2][i,j] + $𝔽[k][i,j] * $m[i,j]
                end
            end
            # $J[k] .= $𝔻[k,1] .* $Du[1] .+ $𝔻[k,2] .* $Du[2] .+ $𝔽[k] .* $m
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


# g_kij, with k the vector component, and ij the grid index
# function wfunc{T}(op::FokkerPlanckOperator{T}, x::T, x_grid::Matrix{T}, g::Vector{Matrix{T}}, w::Matrix{Matrix{T}})
function wfunc_fokker_planck!{T}(g::Vector{Matrix{T}}, w::Matrix{Matrix{T}})
    local g¹::T
    local g²::T

    local t1::Matrix{T}
    local t2::Matrix{T}

    @assert size(w, 1) == length(g)
    @assert size(w, 2) == length(g)
    @assert size(w[1,1], 1) == size(g[1],1)
    @assert size(w[1,1], 2) == size(g[1],2)

    M = size(g[1],1)
    N = size(g[1],2)

    for j in 1:N
        for i in 1:M
            g² = 0
            for k in 1:size(g,1)
                g² += g[k][i,j]^2
                # w[k,k][i,j] = 1
            end
            # g¹ = sqrt(g²)
            for l in 1:size(w,2)
                for k in 1:size(w,1)
                    w[k,l][i,j] = 0
                end
            end
            for k in 1:size(w,1)
                w[k,k][i,j] = g²
            end
            # if g¹ == 0 println("g¹ is zero for [", i, ",", j, "]") end
            for l in 1:size(w,2)
                for k in 1:size(w,1)
                    w[k,l][i,j] -= g[k][i,j] * g[l][i,j]
                    # w[l,k][i,j] -= g[l][i,j] * g[k][i,j] / g²
                    # w[l,k][i,j] /= g¹
                end
            end
        end
    end

    # t1 = zeros(T,M,N)
    # t2 = zeros(T,M,N)
    #
    # for l in 1:length(g)
    #     t1 .+= w[1,l] .* g[l]
    #     t2 .+= w[2,l] .* g[l]
    # end
    #
    # println(t1)
    # println(t2)

    w
end
