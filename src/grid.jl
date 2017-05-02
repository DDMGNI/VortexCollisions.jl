
abstract type Grid{D} end

immutable Grid2d{M,N} <: Grid{2}
    x::Vector
    y::Vector
end

function Grid2d(M, N; x1=0, x2=2π, y1=0, y2=2π)
    x = [x1 + i*(x2-x1)/M for i in 0:M-1]
    y = [y1 + j*(y2-y1)/N for j in 0:N-1]
    Grid2d{M,N}(x, y)
end


@generated function frequencies{M, N, T <: Integer}(r::T, s::T, grid::Grid2d{M,N})
    m = circshift(collect(-div(M-1,2):div(M,2)), -div(M-1,2))
    n = collect(0:div(N,2))

    # println("Creating grid->frequency conversion table for (", M, ",", N, ", grid points:")
    # println("   m = ", m)
    # println("   n = ", n)

    quote
        return $m[r], $n[s]
    end
end
