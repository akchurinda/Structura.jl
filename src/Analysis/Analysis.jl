include("Assembly.jl")
include("Storage.jl")
include("DOFNumberer.jl")
include("DOFPartitioner.jl")
include("TestingSchemes.jl")
include("Algorithms.jl")
include("AnalysisOptions.jl")
include("AnalysisCache.jl")
include("ControlSchemes/LoadControl.jl")
include("ControlSchemes/DisplacementControl.jl")
include("ControlSchemes/WorkControl.jl")

const _NO_LINEAR_LOGGING = LinearVerbosity(LinearSolve.SciMLLogging.None())

function analyze!(model::Model, options::AnalysisOptions, cache::AnalysisCache, steps::Int=1)
    # Check if the number of steps is valid:
    if steps < 1
        throw(ArgumentError("Number of analysis steps must be at least 1."))
    end

    # Extract the analysis options:
    algorithm = options.algorithm
    controlscheme = options.controlscheme
    testingscheme = options.testingscheme

    # Make sure that the algorithm is defined:
    if isnothing(algorithm)
        throw(ArgumentError("Analysis algorithm is not set in AnalysisOptions."))
    end

    # Make sure that the control scheme is defined:
    if isnothing(controlscheme)
        throw(ArgumentError("Analysis control scheme is not set in AnalysisOptions."))
    end

    # Make sure that a testing scheme is provided for the selected algorithm if required:
    if _requires_testing_scheme(algorithm) && isnothing(testingscheme)
        throw(ArgumentError("A testing scheme must be set in AnalysisOptions for nonlinear algorithms."))
    end

    # Extract the committed load factor:
    λ_c = model.λ[]

    # Extract the DOF partitioning:
    fdofs = cache.fdofs

    # Perform the analysis steps:
    converged = true
    for step in 1:steps
        try
            Δλ = _next_step!(model, controlscheme, cache, fdofs)
            λ_t = model.λ[] + Δλ
            
            U, F_res, λ_t = _solve_current_step(controlscheme, model, cache, algorithm, testingscheme, λ_t, fdofs)
            
            _commit_state!(model, cache, U, F_res, λ_t)

            λ_c = λ_t
        catch error
            _revert_state!(model, λ_c)

            @warn "Analysis failed at step $(step)!"

            println("Error message: ", error)

            converged = false
        end
    end

    # Return nothing:
    return converged
end

function _set_trial_state!(model::Model, U::AbstractVector{<:Real})
    # Set the trial state for each node:
    for node in values(model.nodes)
        dofnumbers = _get_dofnumbers(node)
        u_t = SVector{3}(U[dofnumbers[1]], U[dofnumbers[2]], U[dofnumbers[3]])
        _set_trial_state!(node, u_t, node.state.q_t, node.state.r_t)
    end

    # Set the trial state for each element:
    for element in values(model.elements)
        _set_trial_state!(element, _get_u_gbl(element))
    end

    # Set the trial state for the model:
    return nothing
end

function _solve_linear_system!(cache::AnalysisCache, K_red::AbstractMatrix{<:Real}, F_red::AbstractVector{<:Real})
    cache.lscache.A = K_red
    cache.lscache.b = F_red

    solution = solve!(cache.lscache, verbose = _NO_LINEAR_LOGGING)

    return solution.u
end

function _revert_state!(model::Model)
    for node in values(model.nodes)
        _revert_state!(node)
    end
    for element in values(model.elements)
        _revert_state!(element)
    end
    return nothing
end

function _commit_state!(model::Model, cache::AnalysisCache, U::AbstractVector{<:Real}, F_res::AbstractVector{<:Real}, λ_t::Real)
    _commit_state!(model, cache, U, F_res)
    
    model.λ[] = λ_t
    
    return nothing
end

function _commit_state!(model::Model, cache::AnalysisCache, U::AbstractVector{<:Real}, F_res::AbstractVector{<:Real})
    _set_trial_state!(model, U)

    _assemble_Q!(cache.F_int, model)

    for node in values(model.nodes)
        dofnumbers = _get_dofnumbers(node)
        u = SVector{3}(U[dofnumbers[1]], U[dofnumbers[2]], U[dofnumbers[3]])
        q = SVector{3}(cache.F_int[dofnumbers[1]], cache.F_int[dofnumbers[2]], cache.F_int[dofnumbers[3]])
        r = SVector{3}(F_res[dofnumbers[1]], F_res[dofnumbers[2]], F_res[dofnumbers[3]])

        _set_trial_state!(node, u, q, r)
        _commit_state!(node)
    end

    for element in values(model.elements)
        _commit_state!(element)
    end

    return nothing
end

function _revert_state!(model::Model, λ_c::Real)
    _revert_state!(model)

    model.λ[] = λ_c
    
    return nothing
end
