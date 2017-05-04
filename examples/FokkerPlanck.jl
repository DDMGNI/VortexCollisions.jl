
using VortexCollisions

Δt = 1E-6
nt = 10

M = 64
N = 64

Mcut = M
Ncut = N

output = "FokkerPlanck.h5"

gr = Grid2d(M,N)
ft = FourierTransform(gr, mcut=Mcut, ncut=Ncut)
op = FokkerPlanckOperator(gr, ft)


const a = 1
const b = 2

function u_init(x,y)
    exp(a*cos(x-π) + b*cos(y-π))
end


run_simulation(op, nt, Δt, u_init, output)
