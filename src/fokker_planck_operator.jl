
using ParallelDataTransfer


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
        m  = SharedArray{RT}((M,N))

        Dû = [SharedArray{CT}((ℳ,𝒩)), SharedArray{CT}((ℳ,𝒩))]
        Dĥ = [SharedArray{CT}((ℳ,𝒩)), SharedArray{CT}((ℳ,𝒩))]
        m̂  = SharedArray{CT}((ℳ,𝒩))

        𝔽 = [SharedArray{RT}((M,N)), SharedArray{RT}((M,N))]
        𝔻 = Array{SharedArray{RT,2}}(2,2)

        for l in 1:size(𝔻,2)
            for k in 1:size(𝔻,1)
                𝔻[k,l] = SharedArray{RT}((M,N))
            end
        end

        new(grid, ft, D, Δ⁻¹, hfunc, mfunc, wfunc, Du, Dh, m, Dû, Dĥ, m̂, 𝔽, 𝔻)
    end
end

function FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT}(grid::Grid2d{M,N,RT}, ft::FourierTransform{ℳ,𝒩,RT,CT};
                                              hfunc=hfunc_fokker_planck!, mfunc=mfunc_fokker_planck!, wfunc=wfunc_fokker_planck!)
    op = FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT, typeof(hfunc), typeof(mfunc), typeof(wfunc)}(grid, ft, hfunc, mfunc, wfunc)
    nworkers() > 1 ? initialize_workers(op, ft.ℳcut, ft.𝒩cut) : nothing
    return op
end


function initialize_workers{M,N,ℳ,𝒩,RT,CT}(op::FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT}, ℳcut=0, 𝒩cut=0)
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
    local gm  = op.m

    local gDû = op.Dû
    local gm̂  = op.m̂

    local g𝔽  = op.𝔽
    local g𝔻  = op.𝔻

    local gwfunc  = op.wfunc

    local wn::Int = div(N, nworkers())

    @sync for (p,w) in enumerate(workers())
        @spawnat w begin

            global const Dh = gDh
            global const Du = gDu
            global const m  = gm

            global const Dû = gDû
            global const m̂  = gm̂

            global const 𝔽  = g𝔽
            global const 𝔻  = g𝔻

            global const wfunc = gwfunc

            global const j1 = wn*(p-1) + 1
            global const j2 = wn*p

            global const gr = Grid2d(M,N)
            global const ft = FourierTransform(gr, ℳcut=ℳcut, 𝒩cut=𝒩cut)

            # global const lg = [zeros(RT,M,N), zeros(RT,M,N)]
            # global const lw = Array{Array{RT,2}}(2,2)
            # global const lŵ = Array{Array{CT,2}}(2,2)
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


@generated function collision_operator!{M,N,ℳ,𝒩,RT,CT}(op::FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT}, u::Matrix{RT}, divJ::Matrix{RT})

    local ϕ::Matrix{RT}
    local h::Matrix{RT}

    local û::Matrix{CT}
    local ϕ̂::Matrix{CT}
    local ĥ::Matrix{CT}

    local J::Vector{Matrix{RT}}
    local Ĵ::Vector{Matrix{CT}}

    local divĴ::Matrix{CT}


    ϕ  = zeros(RT,M,N)
    h  = zeros(RT,M,N)

    û  = zeros(CT,ℳ,𝒩)
    ϕ̂  = zeros(CT,ℳ,𝒩)
    ĥ  = zeros(CT,ℳ,𝒩)

    J  = [zeros(RT,M,N), zeros(RT,M,N)]
    Ĵ  = [zeros(CT,ℳ,𝒩), zeros(CT,ℳ,𝒩)]

    divĴ = zeros(CT,ℳ,𝒩)


    quote
        @assert size(u) == size(divJ) == (M,N)

        frfft!(op.ft, u, $û)

        apply_operator!(op.Δ⁻¹, $û, $ϕ̂)
        irfft!(op.ft, $ϕ̂, $ϕ)

        op.mfunc(u, op.m)
        frfft!(op.ft, op.m, op.m̂)

        op.hfunc(u, $ϕ, $h)
        frfft!(op.ft, $h, $ĥ)

        apply_operator!(op.D, $û, op.Dû)
        apply_operator!(op.D, $ĥ, op.Dĥ)

        for k in 1:length(op.Dû)
            irfft!(op.ft, op.Dû[k], op.Du[k])
        end

        for k in 1:length(op.Dĥ)
            irfft!(op.ft, op.Dĥ[k], op.Dh[k])
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


function compute_convolution!{M,N,ℳ,𝒩,RT,CT}(op::FokkerPlanckOperator{M,N,ℳ,𝒩,RT,CT})
    if nworkers() > 1
        @sync for w in workers()
            @spawnat w begin
                convolution_kernel_fourier!(gr, ft, j1, j2, wfunc, Dh, Dû, m̂, 𝔽, 𝔻)
                # convolution_kernel_trapezoidal!(gr, j1, j2, wfunc, Dh, Du, m, 𝔽, 𝔻)
            end
        end
    else
        convolution_kernel_fourier!(op.grid, op.ft, 1, N, op.wfunc, op.Dh, op.Dû, op.m̂, op.𝔽, op.𝔻)
        # convolution_kernel_trapezoidal!(gr, 1, N, op.wfunc, op.Dh, op.Du, op.m, op.𝔽, op.𝔻)
    end
end


@generated function convolution_kernel_fourier!{M,N,ℳ,𝒩,RT,CT}(gr::Grid2d{M,N}, ft::FourierTransform{ℳ,𝒩,RT,CT}, j1::Int, j2::Int, wfunc::Function,
                                                                Dh::Vector{SharedArray{RT,2}}, Dû::Vector{SharedArray{CT,2}}, m̂::SharedArray{CT,2},
                                                                𝔽::Vector{SharedArray{RT,2}}, 𝔻::Matrix{SharedArray{RT,2}})

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


@generated function convolution_kernel_trapezoidal!{M,N,RT}(gr::Grid2d{M,N}, j1::Int, j2::Int, wfunc::Function,
                                                            Dh::Vector{SharedArray{RT,2}}, Du::Vector{SharedArray{RT,2}}, m::SharedArray{RT,2},
                                                            𝔽::Vector{SharedArray{RT,2}}, 𝔻::Matrix{SharedArray{RT,2}})

    local g::Vector{Matrix{RT}} = [zeros(RT,M,N), zeros(RT,M,N)]
    local w::Matrix{Matrix{RT}} = Array{Array{RT,2}}(2,2)

    for l in 1:size(w,2)
        for k in 1:size(w,1)
            w[k,l] = zeros(RT,M,N)
        end
    end


    quote
        @assert 1 ≤ j1 ≤ j2 ≤ N

        for j in j1:j2
            for i in 1:M

                for k in 1:length($g)
                    $g[k] .= Dh[k][i,j] .- Dh[k]
                end

                wfunc(gr, i, j, $g, $w)

                for k in 1:length(𝔽)
                    𝔽[k][i,j] = 0
                    for l in 1:size($ŵ,2)
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


function mfunc_fokker_planck!{RT}(u::Matrix{RT}, m::Union{Array{RT,2},SharedArray{RT,2}})
    m .= 1
end


function hfunc_fokker_planck!{RT}(u::Matrix{RT}, ϕ::Matrix{RT}, h::Union{Array{RT,2},SharedArray{RT,2}})
    h .= ϕ
end


@generated function wfunc_fokker_planck!{M,N,RT}(grid::Grid2d{M,N,RT}, i::Int, j::Int, g::Vector{Matrix{RT}}, w::Matrix{Matrix{RT}})
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
