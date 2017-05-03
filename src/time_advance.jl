
@generated function timestep!{RT,CT,M,N,ℳ,𝒩}(op::CollisionOperator{RT,CT,M,N,ℳ,𝒩}, u₀::Matrix{RT}, u₁::Matrix{RT}, Δt::RT)
    local divJ::Matrix{RT} = zeros(RT, M, N)

    quote
        @assert size(u₀) == size(u₁) == size($divJ)
        collision_operator!(op, u₀, $divJ)
        for j in 1:size(u₁,2)
            for i in 1:size(u₁,1)
                u₁[i,j] = u₀[i,j] + Δt * $divJ[i,j]
            end
        end
        u₁
    end
end
