
# to run in parallel set environment variable JULIA_NUM_THREADS to > 1
if haskey(ENV, "JULIA_NUM_THREADS")
    nw = parse(Int, ENV["JULIA_NUM_THREADS"])
    addprocs(nw)
end


using VortexCollisions


const nt = 10000

const M = 128
const N = 128

const ℳcut = 0
const 𝒩cut = 0


r(x,y) = sqrt(x^2 + y^2)

r₀(x,y,ϵ) = 1 + ϵ * cos(2 * Θ(x,y))

Θ(x,y) = atan2(y,x)

function u_init_exp_r(x,y,ϵ=0.4)
    x̅ = (x-π)
    y̅ = (y-π)
    exp(- (r(x̅,y̅) / r₀(x̅,y̅,ϵ))^4)
end


gr = Grid2d(M,N)
ft = FourierTransform(gr, ℳcut=ℳcut, 𝒩cut=𝒩cut)
op = FokkerPlanckOperator(gr, ft)


# test
run_simulation(op, 100, 1E-6, u_init_exp_r,         "FokkerPlanck_exp_r_dt1E-6_nt100.h5")

# production
# run_simulation(op, nt, 1E-6, u_init_exp_r,         "FokkerPlanck_exp_r_dt1E-6_nt"        * string(nt) * ".h5")
