
using Distributed

# to run in parallel set environment variable JULIA_NUM_THREADS to > 1
if haskey(ENV, "JULIA_NUM_THREADS")
    nw = parse(Int, ENV["JULIA_NUM_THREADS"])
    addprocs(nw; exeflags="--project")
end


using VortexCollisions


const nt = 10000

const M = 64
const N = 64

const ℳcut = 0
const 𝒩cut = 0


const a = 1
const b = 2

function u_init_exp_cos(x,y)
    exp(a*cos(x-π) + b*cos(y-π))
end


function u_init_sinx4(x,y)
    sin(x)^4
end


function u_init_sinx4siny4(x,y)
    sin(x)^4 * sin(y)^4
end


function u_init_sinx4siny4etc(x,y)
    sin(x)^4 * sin(y)^4 + 0.2 * sin(2x)^4 * sin(0.5 * y)^4
end



gr = Grid2d(M,N)
ft = FourierTransform(gr, ℳcut=ℳcut, 𝒩cut=𝒩cut)
op = FokkerPlanckOperator(gr, ft)


# test
run_simulation(op, 10, 1E-6, u_init_exp_cos,       "FokkerPlanck_exp_cos_dt1E-6_nt10.h5")
# run_simulation(op, 10, 1E-3, u_init_sinx4,         "FokkerPlanck_sinx4_dt1E-3_nt10.h5")
# run_simulation(op, 10, 1E-3, u_init_sinx4siny4,    "FokkerPlanck_sinx4siny4_dt1E-3_nt10.h5")
# run_simulation(op, 10, 1E-3, u_init_sinx4siny4etc, "FokkerPlanck_sinx4siny4etc_dt1E-3_nt10.h5")

# production
# run_simulation(op, nt, 1E-6, u_init_exp_cos,       "FokkerPlanck_exp_cos_dt1E-6_nt"      * string(nt) * ".h5")
# run_simulation(op, nt, 1e-3, u_init_sinx4,         "FokkerPlanck_sinx4_dt1E-3_nt"        * string(nt) * ".h5")
# run_simulation(op, nt, 1e-3, u_init_sinx4siny4,    "FokkerPlanck_sinx4siny4_dt1E-3_nt"   * string(nt) * ".h5")
# run_simulation(op, nt, 1e-3, u_init_sinx4siny4etc, "FokkerPlanck_sinx4sinyetc_dt1E-3_nt" * string(nt) * ".h5")
