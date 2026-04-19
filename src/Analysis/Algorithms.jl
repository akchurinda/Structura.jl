struct LinearAlgorithm                <: AbstractAlgorithm end
struct NewtonRaphsonAlgorithm         <: AbstractAlgorithm end
struct ModifiedNewtonRaphsonAlgorithm <: AbstractAlgorithm end

_requires_testing_scheme(::LinearAlgorithm               ) = false
_requires_testing_scheme(::NewtonRaphsonAlgorithm        ) = true
_requires_testing_scheme(::ModifiedNewtonRaphsonAlgorithm) = true
