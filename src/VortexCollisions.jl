
module VortexCollisions

    export Grid2d, get_field, frequencies

    include("grid.jl")

    export FourierTransform, get_trans, prfft!, frfft!, irfft!

    include("fourier_transform.jl")

    export apply_operator!, get_gradient, get_inverse_laplacian

    include("fourier_operators.jl")

    export fourier_quadrature

    include("fourier_quadrature.jl")

    export CollisionOperator

    include("collision_operator.jl")

    export FokkerPlanckOperator, collision_operator!

    include("fokker_planck_operator.jl")

    export timestep!

    include("time_integration.jl")

    export run_simulation

    include("simulation.jl")

end
