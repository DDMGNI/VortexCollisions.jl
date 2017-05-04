
"""
Three-stage Runge-Kutta method of Zhao and Wei with parameter C=4
[Math. Meth. Appl. Sci. 2014, 37 1042–1071, figure 1b]
"""
@generated function timestep!{RT,CT,M,N,ℳ,𝒩}(op::CollisionOperator{RT,CT,M,N,ℳ,𝒩}, u₀::Matrix{RT}, u₁::Matrix{RT}, Δt::RT)
    local u::Vector{Matrix{RT}} = [zeros(RT,M,N), zeros(RT,M,N), zeros(RT,M,N)]
    local f::Vector{Matrix{RT}} = [zeros(RT,M,N), zeros(RT,M,N), zeros(RT,M,N)]

    const a::Vector{RT} = [0.0, 0.5, 1.0]
    const b::Vector{RT} = [1/6, 2/3, 1/6]
    const c::Vector{RT} = [0.0, 0.5, 1.0]

    quote
        @assert size(u₀) == size(u₁) == (M,N)

        $u[1] .= u₀
        collision_operator!(op, $u[1], $f[1])

        update_field!($f[1], u₀, $u[2], Δt * $a[2])
        collision_operator!(op, $u[2], $f[2])

        update_field!($f[2], u₀, $u[3], Δt * $a[3])
        collision_operator!(op, $u[3], $f[3])

        update_field!($f[1], u₀, u₁, Δt * $b[1])
        update_field!($f[2], u₁, u₁, Δt * $b[2])
        update_field!($f[3], u₁, u₁, Δt * $b[3])

        u₁
    end
end


function update_field!{RT}(f::Matrix{RT}, u₀::Matrix{RT}, u₁::Matrix{RT}, fac::RT)
    @assert size(f) == size(u₀) == size(u₁)
    @inbounds for j in 1:size(u₁,2)
        for i in 1:size(u₁,1)
            u₁[i,j] = u₀[i,j] + fac * f[i,j]
        end
    end
end