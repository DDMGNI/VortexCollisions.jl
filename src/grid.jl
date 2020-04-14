
abstract type Grid{D} end

struct Grid2d{M,N,T} <: Grid{2}
    x::Vector{T}
    y::Vector{T}
    x1::T
    x2::T
    y1::T
    y2::T
    normalisation::T
end

function Grid2d(M, N; x1=0, x2=2π, y1=0, y2=2π)
    x = [x1 + i*(x2-x1)/M for i in 0:M-1]
    y = [y1 + j*(y2-y1)/N for j in 0:N-1]
    normalisation = (x2 - x1) * (y2 - y1) / M^2 / N^2
    Grid2d{M,N,typeof(normalisation)}(x, y, x1, x2, y1, y2, normalisation)
end


# TODO Add constructor to read grid from HDF5.
# TODO Add function to write grid to HDF5.


function get_field(grid::Grid2d{M,N,T}) where {M,N,T}
    zeros(T, M, N)
end


@generated function frequencies(r::T, s::T, grid::Grid2d{M,N}) where {M, N, T <: Integer}
    m = circshift(collect(-div(M-1,2):div(M,2)), -div(M-1,2))
    n = collect(0:div(N,2))

    # println("Creating grid->frequency conversion table for (", M, ",", N, ", grid points:")
    # println("   m = ", m)
    # println("   n = ", n)

    quote
        return $m[r], $n[s]
    end
end
