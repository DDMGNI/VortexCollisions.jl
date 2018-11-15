
using Distributed

if haskey(ENV, "JULIA_NUM_THREADS")
    nw = parse(Int, ENV["JULIA_NUM_THREADS"])
    addprocs(nw)
else
    @warn "Running in serial mode."
end

const nt = 10
const nx = 16
const ny = 8
const Δt = 1E-1


@everywhere using SharedArrays

@everywhere module CollisionOperator

    using Distributed
    using SharedArrays

    export run_simulation


    struct Grid2d{M,N,T}
        x::Vector{T}
        y::Vector{T}
    end

    function Grid2d(M, N; x1=0, x2=2π, y1=0, y2=2π)
        x = [x1 + i*(x2-x1)/M for i in 0:M-1]
        y = [y1 + j*(y2-y1)/N for j in 0:N-1]
        Grid2d{M,N,eltype(x)}(x, y)
    end


    struct RungeKutta{M,N,T}
        u::Vector{SharedArray{T,2}}
        f::Vector{SharedArray{T,2}}

        Δt::T

        a::Vector{T}
        b::Vector{T}
        c::Vector{T}
    end

    function RungeKutta(M, N, Δt::T) where {T}
        u = [SharedArray{T}((M,N)), SharedArray{T}((M,N)), SharedArray{T}((M,N))]
        f = [SharedArray{T}((M,N)), SharedArray{T}((M,N)), SharedArray{T}((M,N))]

        a = [0.0, 0.5, 1.0]
        b = [1/6, 2/3, 1/6]
        c = [0.0, 0.5, 1.0]

        RungeKutta{M,N,T}(u, f, Δt, a, b, c)
    end


    function run_simulation(nt::Int, nx::Int, ny::Int, Δt::Number, uinit::Function, mfunc::Function)
        grid = Grid2d(nx,ny)
        rk   = RungeKutta(nx, ny, Δt)

        u₀ = zeros(nx,ny)
        u₁ = zeros(nx,ny)

        # compute initial conditions
        for j in 1:size(u₀,2)
            for i in 1:size(u₀,1)
                 u₀[i,j]  = uinit(grid.x[i], grid.y[j])
            end
        end

        # run for nt time steps
        for n in 1:nt
            timestep!(grid, rk, mfunc, u₀, u₁)
            # ...
            # save solution
            # ...
            u₀ .= u₁
        end
    end


    "Runge-Kutta time stepping algorithm solving u'=f(u)"
    function timestep!(grid::Grid2d{M,N,T}, rk::RungeKutta{M,N,T}, mfunc::Function, u₀::Matrix{T}, u₁::Matrix{T}) where {M,N,T}
        @assert size(u₀) == size(u₁) == (M,N)

        # compute intermediate steps
        rk.u[1] .= u₀
        collision_operator!(grid, mfunc, rk.u[1], rk.f[1])

        for i in 2:length(rk.a)
            rk.u[i] .= u₀ .+ rk.Δt .* rk.a[i] .* rk.f[i-1]
            collision_operator!(grid, mfunc, rk.u[i], rk.f[i])
        end

        # compute final solution
        u₁ .= u₀
        for i in eachindex(rk.b, rk.f)
            u₁ .+= rk.Δt .* rk.b[i] .* rk.f[i]
        end
    end


    "Collison operator: computes the vector field f(u) where f is a complicated integro-differential operator"
    function collision_operator!(grid::Grid2d{M,N,T}, m::Function, u::SharedArray{T,2}, f::SharedArray{T,2}) where {M,N,T}
        # ...
        # do something not too expensive
        # ...

        # compute expensive part in parallel
        if nworkers() > 1
            wn = div(N, nworkers())

            # @sync @distributed for j = 1:N
            #     convolution_kernel!(grid, j, j, m, u, f)
            # end

            @sync for (p,w) in enumerate(workers())
                j1 = wn*(p-1) + 1
                j2 = wn*p
                @spawnat w begin
                    convolution_kernel!(grid, j1, j2, m, u, f)
                end
            end

            # @sync for (p,w) in enumerate(workers())
            #     j1 = wn*(p-1) + 1
            #     j2 = wn*p
            #     @async remotecall_wait(convolution_kernel!, w, grid, j1, j2, m, u, f)
            # end

        else
            convolution_kernel!(grid, 1, N, m, u, f)
        end
    end


    "Convolution kernel: most expensive part in vector field computation"
    @generated function convolution_kernel!(grid::Grid2d{M,N,T}, j1::Int, j2::Int, m::Function, u::SharedArray{T,2}, f::SharedArray{T,2}) where {M,N,T}
        # ...
        # create some local arrays
        # ...

        quote
            for j in j1:j2
                for i in 1:M
                    # ...
                    # do some expensive computation involving the function m()
                    # ...
                end
            end
        end
    end

end


# function that is used in the convolution kernel
@everywhere function m_func!(u::Array{T,2}, m::Union{Array{T,2},SharedArray{T,2}}) where {T}
    m .= u
end

# function that sets initial conditions
@everywhere function u_init(x,y)
    sin(x)^4 * sin(y)^4
end

# load CollisionOperator module
using .CollisionOperator

# warmup
run_simulation(nt, nx, ny, Δt, u_init, m_func!)

# timing
@time run_simulation(nt, nx, ny, Δt, u_init, m_func!)
