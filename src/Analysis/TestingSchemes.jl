struct AbsoluteNormUnbalancedForce <: AbstractTestingScheme
    tolerance::Real
    maxiterations::Int
end

struct AbsoluteNormDisplacementIncrement <: AbstractTestingScheme
    tolerance::Real
    maxiterations::Int
end

struct AbsoluteNormEnergyIncrement <: AbstractTestingScheme
    tolerance::Real
    maxiterations::Int
end

_maxiterations(tscheme::AbstractTestingScheme) = tscheme.maxiterations

_compute_test_value(::AbsoluteNormUnbalancedForce      , F_red ::AbstractVector{<:Real}, ΔU_f::AbstractVector{<:Real}) = norm(F_red)
_compute_test_value(::AbsoluteNormDisplacementIncrement, F_red ::AbstractVector{<:Real}, ΔU_f::AbstractVector{<:Real}) = norm(ΔU_f)
_compute_test_value(::AbsoluteNormEnergyIncrement      , F_red ::AbstractVector{<:Real}, ΔU_f::AbstractVector{<:Real}) = abs(dot(ΔU_f, F_red))

_is_converged(tscheme::Nothing, tvalue::Real) = true
_is_converged(tscheme::AbstractTestingScheme, tvalue::Real) = tvalue < tscheme.tolerance