struct DisplacementControl <: AbstractControlScheme
    nid::Int
    dof::Int
    Δu::Real
end

function _next_step!(model::Model, controlscheme::DisplacementControl, cache::AnalysisCache, fdofs::BitVector)
    _assemble_F!(cache.F_ref, model)
    _assemble_K!(cache.K, model)

    reduced_dof = _get_reduced_controlled_dof(model, controlscheme, fdofs)

    _assemble_K_red!(cache.K_red, cache.K, cache.fdof_indices)
    _assemble_F_red!(cache.F_red, cache.F_ref, cache.fdof_indices)

    ΔU_p_f = _solve_linear_system!(cache, cache.K_red, cache.F_red)
    δU_p_c = ΔU_p_f[reduced_dof]

    if iszero(δU_p_c)
        throw(ArgumentError("Displacement Control predictor is singular: controlled DOF has zero reference displacement response."))
    end

    return controlscheme.Δu / δU_p_c
end

function _get_reduced_controlled_dof(model::Model, controlscheme::DisplacementControl, fdofs::BitVector)
    global_dof = _get_dofnumbers(model.nodes[controlscheme.nid])[controlscheme.dof]
    if !fdofs[global_dof]
        throw(ArgumentError("Cannot apply Displacement Control to a constrained DOF."))
    end

    return count(@view fdofs[1:global_dof])
end

function _solve_current_step(
    controlscheme::DisplacementControl,
    model::Model,
    cache::AnalysisCache,
    algorithm::LinearAlgorithm,
    testingscheme::Union{Nothing,AbstractTestingScheme},
    λ_t::Real,
    fdofs::BitVector)

    _assemble_U!(cache.U, model)
    _assemble_F!(cache.F_ref, model)
    _assemble_K!(cache.K, model)

    reduced_dof = _get_reduced_controlled_dof(model, controlscheme, fdofs)

    _assemble_K_red!(cache.K_red, cache.K, cache.fdof_indices)
    _assemble_F_red!(cache.F_red, cache.F_ref, cache.fdof_indices)

    # Solve for reference displacements
    ΔU_p_f = _solve_linear_system!(cache, cache.K_red, cache.F_red)

    # Apply predictor increment implied by _next_step!:
    Δλ = λ_t - model.λ[]

    # Scale and update global displacements:
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
    controlscheme::DisplacementControl,
    model::Model,
    cache::AnalysisCache,
    algorithm::Union{NewtonRaphsonAlgorithm,ModifiedNewtonRaphsonAlgorithm},
    testingscheme::Union{Nothing,AbstractTestingScheme},
    λ_t::Real,
    fdofs::BitVector)

    _assemble_U!(cache.U, model)
    _assemble_F!(cache.F_ref, model)

    reduced_dof = _get_reduced_controlled_dof(model, controlscheme, fdofs)

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
        δU_p_f = _solve_linear_system!(cache, cache.K_red, cache.F_red)
        
        # Store securely in the cache to avoid LinearSolve.jl memory overwrites
        @inbounds for (i_local, i_global) in enumerate(cache.fdof_indices)
            cache.δU_p[i_global] = δU_p_f[i_local]
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
        δU_p_c = cache.δU_p[cache.fdof_indices[reduced_dof]]
        
        if i == 1
            # Predictor λ increment was computed by _next_step!.
            δλ = λ_t - model.λ[]
        else
            δU_r_c = cache.δU_r[cache.fdof_indices[reduced_dof]]
            δλ = -δU_r_c / δU_p_c
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

        testvalue = _compute_test_value(testingscheme, cache.F_red, δU_p_f)

        if _is_converged(testingscheme, testvalue)
            return cache.U, cache.F_res, λ_t # Return the updated λ_t!
        end
    end

    throw(ErrorException("Nonlinear iterations did not converge in $(maxiterations) iterations."))
end