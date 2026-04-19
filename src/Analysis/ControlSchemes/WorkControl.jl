using LinearAlgebra

struct WorkControl <: AbstractControlScheme
    ΔW::Real
end

function _next_step!(model::Model, controlscheme::WorkControl, cache::AnalysisCache, fdofs::BitVector)
    _assemble_F!(cache.F_ref, model)
    _assemble_K!(cache.K, model)

    _assemble_K_red!(cache.K_red, cache.K, cache.fdof_indices)
    _assemble_F_red!(cache.F_red, cache.F_ref, cache.fdof_indices)

    ΔU_p_f = _solve_linear_system!(cache, cache.K_red, cache.F_red)
    
    # Calculate reference work
    W_ref = dot(cache.F_red, ΔU_p_f)

    if iszero(W_ref)
        throw(ArgumentError("Work Control predictor is singular: reference work is zero."))
    end

    # Predictor: use sign(W_ref) to automatically flip the loading direction past limit points
    return sign(W_ref) * sqrt(abs(controlscheme.ΔW / W_ref))
end

function _solve_current_step(
    controlscheme::WorkControl,
    model::Model,
    cache::AnalysisCache,
    algorithm::LinearAlgorithm,
    testingscheme::Union{Nothing,AbstractTestingScheme},
    λ_t::Real,
    fdofs::BitVector)

    _assemble_U!(cache.U, model)
    _assemble_F!(cache.F_ref, model)
    _assemble_K!(cache.K, model)

    _assemble_K_red!(cache.K_red, cache.K, cache.fdof_indices)
    _assemble_F_red!(cache.F_red, cache.F_ref, cache.fdof_indices)

    # Solve for reference displacements
    ΔU_p_f = _solve_linear_system!(cache, cache.K_red, cache.F_red)

    # Apply predictor increment implied by _next_step!:
    Δλ = λ_t - model.λ[]

    # Scale and update global displacements
    @inbounds for (i_local, i_global) in enumerate(cache.fdof_indices)
        cache.U[i_global] += Δλ * ΔU_p_f[i_local]
    end

    _set_trial_state!(model, cache.U)
    _assemble_Q!(cache.F_int, model)

    @inbounds cache.F_ext .= λ_t .* cache.F_ref
    @inbounds cache.F_res .= cache.F_ext .- cache.F_int

    # Return the updated λ_t!
    return cache.U, cache.F_res, λ_t 
end

function _solve_current_step(
    controlscheme::WorkControl,
    model::Model,
    cache::AnalysisCache,
    algorithm::Union{NewtonRaphsonAlgorithm,ModifiedNewtonRaphsonAlgorithm},
    testingscheme::Union{Nothing,AbstractTestingScheme},
    λ_t::Real,
    fdofs::BitVector)

    _assemble_U!(cache.U, model)
    _assemble_F!(cache.F_ref, model)

    maxiterations = _maxiterations(testingscheme)
    recompute = algorithm isa NewtonRaphsonAlgorithm
    testvalue = zero(eltype(cache.U))

    _set_trial_state!(model, cache.U)
    _assemble_Q!(cache.F_int, model)
    @inbounds cache.F_ext .= λ_t .* cache.F_ref
    @inbounds cache.F_res .= cache.F_ext .- cache.F_int

    _assemble_K!(cache.K, model)
    _assemble_K_red!(cache.K_red, cache.K, cache.fdof_indices)

    for i in 1:maxiterations
        if i > 1 && (recompute || i == 1)
            _assemble_K!(cache.K, model)
            _assemble_K_red!(cache.K_red, cache.K, cache.fdof_indices)
        end

        # --- Sub-solve 1: Tangent displacement due to Reference Loads ---
        _assemble_F_red!(cache.F_red, cache.F_ref, cache.fdof_indices)
        
        # Only re-solve for the reference load if the stiffness matrix changed
        if recompute || i == 1
            δU_p_f = _solve_linear_system!(cache, cache.K_red, cache.F_red)
            @inbounds for (i_local, i_global) in enumerate(cache.fdof_indices)
                cache.δU_p[i_global] = δU_p_f[i_local]
            end
        end

        # --- Sub-solve 2: Tangent displacement due to Unbalanced Forces ---
        if i == 1
            # Zero the residual displacements on first iteration
            for j_idx in cache.fdof_indices
                cache.δU_r[j_idx] = zero(eltype(cache.U))
            end
        else
            _assemble_F_red!(cache.F_red, cache.F_res, cache.fdof_indices)
            δU_r_f = _solve_linear_system!(cache, cache.K_red, cache.F_red)
            @inbounds for (i_local, i_global) in enumerate(cache.fdof_indices)
                cache.δU_r[i_global] = δU_r_f[i_local]
            end
        end

        # --- Enforce Constraint & Calculate δλ ---
        if i == 1
            # Predictor λ increment was computed by _next_step!.
            δλ = λ_t - model.λ[]
        else
            # Corrector step: make the iterative work done exactly zero
            _assemble_F_red!(cache.F_red, cache.F_ref, cache.fdof_indices)
            
            # Calculate work contributions
            W_ref_local = zero(eltype(cache.U))
            W_res_local = zero(eltype(cache.U))
            
            @inbounds for (i_local, i_global) in enumerate(cache.fdof_indices)
                W_ref_local += cache.F_red[i_local] * cache.δU_p[i_global]
                W_res_local += cache.F_red[i_local] * cache.δU_r[i_global]
            end
            
            δλ = -W_res_local / W_ref_local
            λ_t += δλ
        end

        # --- Calculate Total Displacement Increment & Update ---
        for (i_local, i_global) in enumerate(cache.fdof_indices)
            cache.U[i_global] += cache.δU_r[i_global] + δλ * cache.δU_p[i_global]
        end

        # --- Re-evaluate State ---
        _set_trial_state!(model, cache.U)
        _assemble_Q!(cache.F_int, model)
        @inbounds cache.F_ext .= λ_t .* cache.F_ref
        @inbounds cache.F_res .= cache.F_ext .- cache.F_int
        _assemble_F_red!(cache.F_red, cache.F_res, cache.fdof_indices)

        # Compute reduced-system displacement increment for convergence check
        @inbounds for (i_local, i_global) in enumerate(cache.fdof_indices)
            cache.U_f[i_local] = cache.δU_r[i_global] + δλ * cache.δU_p[i_global]
        end
        testvalue = _compute_test_value(testingscheme, cache.F_red, cache.U_f)

        if _is_converged(testingscheme, testvalue)
            return cache.U, cache.F_res, λ_t # Return the updated λ_t!
        end
    end

    throw(ErrorException("Nonlinear iterations did not converge in $(maxiterations) iterations."))
end
