struct LoadControl <: AbstractControlScheme
    Δλ::Real
end

function _next_step!(model::Model, controlscheme::LoadControl, ::AnalysisCache, ::BitVector)
    return controlscheme.Δλ
end

function _solve_current_step(
    ::LoadControl,
    model::Model,
    cache::AnalysisCache,
    algorithm::LinearAlgorithm,
    testingscheme::Union{Nothing,AbstractTestingScheme},
    λ_t::Real,
    fdofs::BitVector)

    _assemble_U!(cache.U, model)

    _assemble_F!(cache.F_ref, model)
    @inbounds cache.F_ext .= λ_t .* cache.F_ref

    _set_trial_state!(model, cache.U)
    _assemble_Q!(cache.F_int, model)
    @inbounds cache.F_res .= cache.F_ext .- cache.F_int

    _assemble_K!(cache.K, model)
    _assemble_K_red!(cache.K_red, cache.K, cache.fdof_indices)
    _assemble_F_red!(cache.F_red, cache.F_res, cache.fdof_indices)

    ΔU_f = _solve_linear_system!(cache, cache.K_red, cache.F_red)

    @inbounds for (i_local, i_global) in enumerate(cache.fdof_indices)
        cache.U[i_global] += ΔU_f[i_local]
    end

    _set_trial_state!(model, cache.U)
    _assemble_Q!(cache.F_int, model)
    @inbounds cache.F_res .= cache.F_ext .- cache.F_int

    return cache.U, cache.F_res, λ_t
end

function _solve_current_step(
    ::LoadControl,
    model::Model,
    cache::AnalysisCache,
    algorithm::Union{NewtonRaphsonAlgorithm,ModifiedNewtonRaphsonAlgorithm},
    testingscheme::Union{Nothing,AbstractTestingScheme},
    λ_t::Real,
    fdofs::BitVector)

    _assemble_U!(cache.U, model)

    _assemble_F!(cache.F_ref, model)
    @inbounds cache.F_ext .= λ_t .* cache.F_ref

    maxiterations = _maxiterations(testingscheme)
    recompute = algorithm isa NewtonRaphsonAlgorithm

    testvalue = zero(eltype(cache.U))

    _set_trial_state!(model, cache.U)
    _assemble_Q!(cache.F_int, model)
    @inbounds cache.F_res .= cache.F_ext .- cache.F_int
    _assemble_F_red!(cache.F_red, cache.F_res, cache.fdof_indices)

    _assemble_K!(cache.K, model)
    _assemble_K_red!(cache.K_red, cache.K, cache.fdof_indices)

    for i in 1:maxiterations
        if i > 1 && (recompute || i == 1)
            _assemble_K!(cache.K, model)
            _assemble_K_red!(cache.K_red, cache.K, cache.fdof_indices)
        end

        ΔU_f = _solve_linear_system!(cache, cache.K_red, cache.F_red)
        
        @inbounds for (i_local, i_global) in enumerate(cache.fdof_indices)
            cache.U[i_global] += ΔU_f[i_local]
        end

        _set_trial_state!(model, cache.U)
        _assemble_Q!(cache.F_int, model)
        @inbounds cache.F_res .= cache.F_ext .- cache.F_int
        _assemble_F_red!(cache.F_red, cache.F_res, cache.fdof_indices)

        testvalue = _compute_test_value(testingscheme, cache.F_red, ΔU_f)

        if _is_converged(testingscheme, testvalue)
            @inbounds cache.F_res .= cache.F_ext .- cache.F_int
            # Add λ_t to the return tuple
            return cache.U, cache.F_res, λ_t
        end
    end

    throw(ErrorException("Nonlinear maxiterations did not converge in $(maxiterations) maxiterations."))
end