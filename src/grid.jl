
abstract type Grid{D} end

immutable Grid2d{M,N} <: Grid{2}
end

function Grid2d(M,N)
    Grid2d{M,N}()
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
