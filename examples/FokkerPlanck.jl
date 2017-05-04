
using VortexCollisions

Δt = 1E-6
nt = 10

M = 64
N = 64

ℳcut = M
𝒩cut = N

output = "FokkerPlanck.h5"


const a = 1
const b = 2

function u_init(x,y)
    exp(a*cos(x-π) + b*cos(y-π))
end


gr = Grid2d(M,N)
ft = FourierTransform(gr, ℳcut=ℳcut, 𝒩cut=𝒩cut)
op = FokkerPlanckOperator(gr, ft)


run_simulation(op, nt, Δt, u_init, output)
