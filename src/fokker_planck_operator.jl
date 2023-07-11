

struct FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT,HT,MT,WT} <: CollisionOperator{M,N,ℳ,𝒩,RT,CT}
    grid::Grid2d{M,N,RT}
    ft::FourierTransform{ℳ,𝒩,RT,CT}
    D::Vector{Matrix{CT}}
    Δ⁻¹::Matrix{CT}

    hfunc::HT
    mfunc::MT
    wfunc::WT

    Du::Vector{SharedArray{RT,2}}
    Dh::Vector{SharedArray{RT,2}}
    h::SharedArray{RT,2}
    m::SharedArray{RT,2}

    Dû::Vector{SharedArray{CT,2}}
    Dĥ::Vector{SharedArray{CT,2}}
    m̂::SharedArray{CT,2}

    𝔽::Vector{SharedArray{RT,2}}
    𝔻::Matrix{SharedArray{RT,2}}


    function FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT,HT,MT,WT}(grid::Grid2d{M,N,RT}, ft::FourierTransform{ℳ,𝒩,RT,CT},
                                                           hfunc::HT, mfunc::MT, wfunc::WT) where {M,N,ℳ,𝒩,RT,CT,HT,MT,WT}
        D   = get_gradient(ft)
        Δ⁻¹ = get_inverse_laplacian(ft, sign=-1)

        Du = [SharedArray{RT}((M,N)), SharedArray{RT}((M,N))]
        Dh = [SharedArray{RT}((M,N)), SharedArray{RT}((M,N))]
        h  = SharedArray{RT}((M,N))
        m  = SharedArray{RT}((M,N))

        Dû = [SharedArray{CT}((ℳ,𝒩)), SharedArray{CT}((ℳ,𝒩))]
        Dĥ = [SharedArray{CT}((ℳ,𝒩)), SharedArray{CT}((ℳ,𝒩))]
        m̂  = SharedArray{CT}((ℳ,𝒩))

        𝔽 = [SharedArray{RT}((M,N)), SharedArray{RT}((M,N))]
        𝔻 = [SharedArray{RT}((M,N)) for k ∈ 1:2, l ∈ 1:2]

        new(grid, ft, D, Δ⁻¹, hfunc, mfunc, wfunc, Du, Dh, h, m, Dû, Dĥ, m̂, 𝔽, 𝔻)
    end
end

function FokkerPlanckOperator(grid::Grid2d{M,N,RT}, ft::FourierTransform{ℳ,𝒩,RT,CT};
                              hfunc=hfunc_fokker_planck!, mfunc=mfunc_fokker_planck!, wfunc=wfunc_fokker_planck!) where {M,N,ℳ,𝒩,RT,CT}
    op = FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT, typeof(hfunc), typeof(mfunc), typeof(wfunc)}(grid, ft, hfunc, mfunc, wfunc)
    nworkers() > 1 ? initialize_workers(op, ft.ℳcut, ft.𝒩cut) : nothing
    return op
end


function initialize_workers(op::FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT}, ℳcut=0, 𝒩cut=0) where {M,N,ℳ,𝒩,RT,CT}
    # # add parallel workers
    # if haskey(ENV, "JULIA_NUM_THREADS")
    #     nw = parse(Int, ENV["JULIA_NUM_THREADS"])
    # else
    #     nw = 1
    # end
    #
    # addprocs(nw)

    local gDh = op.Dh
    local gDu = op.Du
    local gh  = op.h
    local gm  = op.m

    local gDû = op.Dû
    local gm̂  = op.m̂

    local g𝔽  = op.𝔽
    local g𝔻  = op.𝔻

    local gwfunc  = op.wfunc

    local wn::Int = div(N, nworkers())

    # @eval @everywhere Dh = $gDh
    # @eval @everywhere Du = $gDu
    # @eval @everywhere m  = $gm

    @sync for (p,w) in enumerate(workers())
        @spawnat w begin

            global Dh = gDh
            global Du = gDu
            global h  = gh
            global m  = gm

            global Dû = gDû
            global m̂  = gm̂

            global 𝔽  = g𝔽
            global 𝔻  = g𝔻

            global wfunc = gwfunc

            global j1 = wn*(p-1) + 1
            global j2 = wn*p

            global gr = Grid2d(M,N)
            global ft = FourierTransform(gr, ℳcut=ℳcut, 𝒩cut=𝒩cut)

            # global lg = [zeros(RT,M,N), zeros(RT,M,N)]
            # global lw = Array{Array{RT,2}}(2,2)
            # global lŵ = Array{Array{CT,2}}(2,2)
            #
            # for l in 1:size(lw,2)
            #     for k in 1:size(lw,1)
            #         lw[k,l] = zeros(RT,M,N)
            #     end
            # end
            #
            # for l in 1:size(lŵ,2)
            #     for k in 1:size(lŵ,1)
            #         lŵ[k,l] = zeros(CT,ℳ,𝒩)
            #     end
            # end

        end
    end
end


@generated function collision_operator!(op::FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT,HT,MT,WT}, u::Matrix{RT}, divJ::Matrix{RT}) where {M,N,ℳ,𝒩,RT,CT,HT,MT,WT}

    local ϕ::Matrix{RT}
    local h::Matrix{RT}
    local m::Matrix{RT}

    local Du::Matrix{RT}
    local Dh::Matrix{RT}

    local û::Matrix{CT}
    local ϕ̂::Matrix{CT}
    local ĥ::Matrix{CT}
    local m̂::Matrix{CT}

    local Dû::Matrix{CT}
    local Dĥ::Matrix{CT}

    local J::Vector{Matrix{RT}}
    local Ĵ::Vector{Matrix{CT}}

    local divĴ::Matrix{CT}


    ϕ  = zeros(RT,M,N)
    h  = zeros(RT,M,N)
    m  = zeros(RT,M,N)

    Du = zeros(RT,M,N)
    Dh = zeros(RT,M,N)

    û  = zeros(CT,ℳ,𝒩)
    ϕ̂  = zeros(CT,ℳ,𝒩)
    ĥ  = zeros(CT,ℳ,𝒩)
    m̂  = zeros(CT,ℳ,𝒩)

    Dû = zeros(CT,ℳ,𝒩)
    Dĥ = zeros(CT,ℳ,𝒩)

    J  = [zeros(RT,M,N), zeros(RT,M,N)]
    Ĵ  = [zeros(CT,ℳ,𝒩), zeros(CT,ℳ,𝒩)]

    divĴ = zeros(CT,ℳ,𝒩)


    quote
        @assert size(u) == size(divJ) == (M,N)

        frfft!(op.ft, u, $û)

        apply_operator!(op.Δ⁻¹, $û, $ϕ̂)
        irfft!(op.ft, $ϕ̂, $ϕ)

        op.mfunc(u, $m, op.grid)
        frfft!(op.ft, $m, $m̂)

        op.m .= $m
        op.m̂ .= $m̂

        op.hfunc(u, $ϕ, $h, op.grid)
        frfft!(op.ft, $h, $ĥ)

        op.h .= $h
        
        apply_operator!(op.D, $û, op.Dû)
        apply_operator!(op.D, $ĥ, op.Dĥ)

        # FFT does not seem to work on Shared Array!
        #
        # for k in 1:length(op.Du)
        #     irfft!(op.ft, op.Dû[k], op.Du[k])
        # end
        #
        # for k in 1:length(op.Dh)
        #     irfft!(op.ft, op.Dĥ[k], op.Dh[k])
        # end

        for k in 1:length(op.Du)
            $Dû .= op.Dû[k]
            irfft!(op.ft, $Dû, $Du)
            op.Du[k] .= $Du
        end

        for k in 1:length(op.Dh)
            $Dĥ .= op.Dĥ[k]
            irfft!(op.ft, $Dĥ, $Dh)
            op.Dh[k] .= $Dh
        end


        compute_convolution!(op)

        for k in 1:length($J)
            @inbounds for j in 1:size($J[k],2)
                for i in 1:size($J[k],1)
                    $J[k][i,j] = op.𝔽[k][i,j] * op.m[i,j]
                end
            end

            for l in 1:length($J)
                @inbounds for j in 1:size($J[k],2)
                    for i in 1:size($J[k],1)
                        $J[k][i,j] += op.𝔻[k,l][i,j] * op.Du[l][i,j]
                    end
                end
            end

            prfft!(op.ft, $J[k], $Ĵ[k])
        end

        apply_operator!(op.D, $Ĵ, $divĴ)

        irfft!(op.ft, $divĴ, divJ)

        divJ
    end
end


function compute_convolution!(op::FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT}) where {M,N,ℳ,𝒩,RT,CT}
    if nworkers() > 1
        # @sync @parallel for j = 1:N
        #     # convolution_kernel_fourier!(gr, ft, j, j, wfunc, Dh, Dû, m̂, 𝔽, 𝔻)
        #     convolution_kernel_trapezoidal!(gr, j, j, wfunc, Dh, Du, m, 𝔽, 𝔻)
        # end

        @sync for w in workers()
            @spawnat w begin
                # convolution_kernel_fourier!(gr, ft, j1, j2, wfunc, Dh, Dû, m̂, 𝔽, 𝔻)
                convolution_kernel_trapezoidal!(gr, j1, j2, wfunc, Dh, Du, m, 𝔽, 𝔻)
            end
        end

        # @sync for w in workers()
        #     # @async remotecall_wait(convolution_kernel_fourier!, w, gr, ft, j1, j2, wfunc, Dh, Dû, m̂, 𝔽, 𝔻)
        #     @async remotecall_wait(convolution_kernel_trapezoidal!, w, gr, j1, j2, wfunc, Dh, Du, m, 𝔽, 𝔻)
        # end

    else
        # convolution_kernel_fourier!(op.grid, op.ft, 1, N, op.wfunc, op.Dh, op.Dû, op.m̂, op.𝔽, op.𝔻)
        convolution_kernel_trapezoidal!(op.grid, 1, N, op.wfunc, op.Dh, op.Du, op.m, op.𝔽, op.𝔻)
    end
end


@generated function convolution_kernel_fourier!(gr::Grid2d{M,N}, ft::FourierTransform{ℳ,𝒩,RT,CT}, j1::Int, j2::Int, wfunc::Function,
                                                Dh::Vector{SharedArray{RT,2}}, Dû::Vector{SharedArray{CT,2}}, m̂::SharedArray{CT,2},
                                                𝔽::Vector{SharedArray{RT,2}}, 𝔻::Matrix{SharedArray{RT,2}}) where {M,N,ℳ,𝒩,RT,CT}

    local g::Vector{Matrix{RT}} = [zeros(RT,M,N), zeros(RT,M,N)]
    local w::Matrix{Matrix{RT}} = Array{Array{RT,2}}(2,2)
    local ŵ::Matrix{Matrix{CT}} = Array{Array{CT,2}}(2,2)

    for l in 1:size(w,2)
        for k in 1:size(w,1)
            w[k,l] = zeros(RT,M,N)
        end
    end

    for l in 1:size(ŵ,2)
        for k in 1:size(ŵ,1)
            ŵ[k,l] = zeros(CT,ℳ,𝒩)
        end
    end


    quote
        @assert 1 ≤ j1 ≤ j2 ≤ N
        @assert length($g) == length(Dh) == length(Dû)
        @assert length($g) == length(𝔽) == size(𝔻,1) == size(𝔻,2)

        for j in j1:j2
            for i in 1:M

                for k in 1:length($g)
                    $g[k] .= Dh[k][i,j] .- Dh[k]
                end

                wfunc(gr, i, j, $g, $w)

                for l in 1:size($w,2)
                    for k in 1:size($w,1)
                        prfft!(ft, $w[k,l], $ŵ[k,l])
                    end
                end

                for k in 1:length(𝔽)
                    𝔽[k][i,j] = 0
                    for l in 1:length(𝔽)
                        𝔽[k][i,j] -= real(fourier_quadrature($ŵ[k,l], Dû[l], ft.μ, gr))
                    end
                end

                for l in 1:size(𝔻,2)
                    for k in 1:size(𝔻,1)
                        𝔻[k,l][i,j] = real(fourier_quadrature($ŵ[k,l], m̂, ft.μ, gr))
                    end
                end

            end
        end
    end
end


@generated function convolution_kernel_trapezoidal!(gr::Grid2d{M,N}, j1::Int, j2::Int, wfunc::Function,
                                                    Dh::Vector{SharedArray{RT,2}}, Du::Vector{SharedArray{RT,2}}, m::SharedArray{RT,2},
                                                    𝔽::Vector{SharedArray{RT,2}}, 𝔻::Matrix{SharedArray{RT,2}}) where {M,N,RT}

    local g::Vector{Matrix{RT}} = [zeros(RT,M,N), zeros(RT,M,N)]
    local w::Matrix{Matrix{RT}} = [zeros(RT,M,N) for k ∈ 1:2, l ∈ 1:2]


    quote
        @assert 1 ≤ j1 ≤ j2 ≤ N
        @assert length($g) == length(Dh) == length(Du)
        @assert length($g) == length(𝔽) == size(𝔻,1) == size(𝔻,2)

        for j in j1:j2
            for i in 1:M

                for k in 1:length($g)
                    $g[k] .= Dh[k][i,j] .- Dh[k]
                end

                wfunc(gr, i, j, $g, $w)

                for k in 1:length(𝔽)
                    𝔽[k][i,j] = 0
                    for l in 1:length(𝔽)
                        𝔽[k][i,j] -= trapezoidal_quadrature($w[k,l], Du[l], gr)
                    end
                end

                for l in 1:size(𝔻,2)
                    for k in 1:size(𝔻,1)
                        𝔻[k,l][i,j] = trapezoidal_quadrature($w[k,l], m, gr)
                    end
                end

            end
        end
    end
end


function mfunc_fokker_planck!(u::Matrix{RT}, m::Union{Array{RT,2},SharedArray{RT,2}}, grid::Grid2d) where {RT}
    m .= 1
end


function hfunc_fokker_planck!(u::Matrix{RT}, ϕ::Matrix{RT}, h::Union{Array{RT,2},SharedArray{RT,2}}, grid::Grid2d) where {RT}
    h .= ϕ
end


@generated function wfunc_fokker_planck!(grid::Grid2d{M,N,RT}, i::Int, j::Int, g::Vector{Matrix{RT}}, w::Matrix{Matrix{RT}}) where {M,N,RT}
    # local g¹::Matrix{RT} = zeros(RT,M,N)
    local g²::Matrix{RT} = zeros(RT,M,N)

    quote
        @assert size(w,1) == size(w,2) == length(g)
        @assert size(w[1,1],1) == size(g[1],1) == M
        @assert size(w[1,1],2) == size(g[1],2) == N

        for l in 1:size(w,2)
            for k in 1:size(w,1)
                w[k,l] .= 0
            end
        end

        $g² .= 0

        for k in 1:size(g,1)
            @inbounds for j in 1:size($g²,2)
                for i in 1:size($g²,1)
                    $g²[i,j] += g[k][i,j]^2
                end
            end
        end

        # $g¹ .= sqrt.($g²)

        for k in 1:size(w,1)
            @inbounds for j in 1:size(w[k,k],2)
                for i in 1:size(w[k,k],1)
                    w[k,k][i,j] = $g²[i,j]
                end
            end
        end

        for l in 1:size(w,2)
            for k in 1:size(w,1)
                @inbounds for j in 1:size(w[k,l],2)
                    for i in 1:size(w[k,l],1)
                        w[k,l][i,j] -= g[k][i,j] * g[l][i,j]
                    end
                end
            end
        end

        w
    end
end
