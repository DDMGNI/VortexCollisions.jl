
if haskey(ENV, "JULIA_NUM_THREADS")
    nw = parse(Int, ENV["JULIA_NUM_THREADS"])
    addprocs(nw)
end


using VortexCollisions


const Δt = 1E-6
const nt = 2

const M = 64
const N = 64

const ℳcut = M
const 𝒩cut = N

const output = "FokkerPlanck.h5"


const a = 1
const b = 2

function u_init(x,y)
    exp(a*cos(x-π) + b*cos(y-π))
end


gr = Grid2d(M,N)
ft = FourierTransform(gr, ℳcut=ℳcut, 𝒩cut=𝒩cut)
op = FokkerPlanckOperator(gr, ft)

run_simulation(op, nt, Δt, u_init, output)
