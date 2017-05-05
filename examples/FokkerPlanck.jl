
# to run in parallel set environment variable JULIA_NUM_THREADS to > 1
if haskey(ENV, "JULIA_NUM_THREADS")
    nw = parse(Int, ENV["JULIA_NUM_THREADS"])
    addprocs(nw)
end


using VortexCollisions


const Δt = 1E-6
const nt = 10000

const M = 64
const N = 64

const ℳcut = M
const 𝒩cut = N


const a = 1
const b = 2

function u_init_exp_cos(x,y)
    exp(a*cos(x-π) + b*cos(y-π))
end


function u_init_sinx4(x,y)
    sin(x)^4
end


function u_init_sin4x(x,y)
    sin(4x)
end


gr = Grid2d(M,N)
ft = FourierTransform(gr, ℳcut=ℳcut, 𝒩cut=𝒩cut)
op = FokkerPlanckOperator(gr, ft)


# test
# run_simulation(op, 10, Δt, u_init_exp_cos, "FokkerPlanck_exp_cos_dt1E-6_nt10.h5")
# run_simulation(op, 10, Δt, u_init_sinx4,   "FokkerPlanck_sinx4_dt1E-6_nt10.h5")
# run_simulation(op, 10, Δt, u_init_sin4x,   "FokkerPlanck_sin4x_dt1E-6_nt10.h5")

# production
# run_simulation(op, nt, Δt, u_init_exp_cos, "FokkerPlanck_exp_cos_dt1E-6_nt" * string(nt) * ".h5")
# run_simulation(op, nt, Δt, u_init_sinx4,   "FokkerPlanck_sinx4_dt1E-6_nt" * string(nt) * ".h5")
# run_simulation(op, nt, Δt, u_init_sin4x,   "FokkerPlanck_sin4x_dt1E-6_nt" * string(nt) * ".h5")
