
using Distributed
using SharedArrays

# to run in parallel set environment variable JULIA_NUM_THREADS to > 1
if haskey(ENV, "JULIA_NUM_THREADS")
    nw = parse(Int, ENV["JULIA_NUM_THREADS"])
    addprocs(nw; exeflags = "--project")
end

using VortexCollisions

const nt = 25000
const Δt = 2E-6

const M = 64
const N = 64

const ℳcut = 0
const 𝒩cut = 0

const a = 1
const b = 2

function u_init(x, y)
    1 + sin(x)^4 * sin(y)^4
end

function mfunc_landau!(u::Matrix{RT}, m::Union{Array{RT, 2}, SharedArray{RT, 2}}, grid::Grid2d) where {RT}
    m .= u
end

function hfunc_landau!(u::Matrix{RT}, ϕ::Matrix{RT},
        h::Union{Array{RT, 2}, SharedArray{RT, 2}}, grid::Grid2d) where {RT}
    for i in axes(h, 1)
        for j in axes(h, 2)
            h[i, j] = (grid.x[i]^2 + grid.y[j]^2) / 2
        end
    end
end

gr = Grid2d(M, N; x1 = -π, x2 = +π, y1 = -π, y2 = +π)
ft = FourierTransform(gr, ℳcut = ℳcut, 𝒩cut = 𝒩cut)
op = FokkerPlanckOperator(gr, ft; hfunc = hfunc_landau!, mfunc = mfunc_landau!)#, wfunc=wfunc_landau!)

# test
# run_simulation(op, 10, Δt, u_init, "Landau_dt1E-3_t10.h5")

# production
# run_simulation(op, nt, Δt, u_init, "Landau_dt1E-6_nt" * string(nt) * ".h5")
run_simulation(op, nt, Δt, u_init, "Landau_dt2E-6_nt" * string(nt) * ".h5")
