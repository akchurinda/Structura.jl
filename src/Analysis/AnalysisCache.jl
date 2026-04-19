struct AnalysisCache{L}
    fdofs::BitVector
    sdofs::BitVector
    fdof_indices::Vector{Int}
    sdof_indices::Vector{Int}
    F_ref::AbstractVector{<:Real}
    F_ext::AbstractVector{<:Real}
    F_int::AbstractVector{<:Real}
    F_res::AbstractVector{<:Real}
    U    ::AbstractVector{<:Real}
    δU_p ::AbstractVector{<:Real} # Iterative displacement due to reference load
    δU_r ::AbstractVector{<:Real} # Iterative displacement due to residual load
    K    ::AbstractMatrix{<:Real}
    K_red::AbstractMatrix{<:Real}
    F_red::AbstractVector{<:Real}
    U_red::AbstractVector{<:Real}
    lscache::L
end

function prepare(model::Model, options::AnalysisOptions)
    # Determine the type of the model:
    T = _get_T(model)

    # Compute the total number of DOFs:
    num_dofs  = 3 * length(model.nodes)

    # Number the DOFs according to the provided DOF numberer:
    numberdofs!(model, options.dofnumberer)

    # Partition the DOFs into free and supported sets:
    fdofs, sdofs = partitiondofs(model)

    # Compute the number of free and supported DOFs:
    num_fdofs = count(fdofs)
    num_sdofs = count(sdofs)

    # Get the indices of the free and supported DOFs:
    fdof_indices = findall(fdofs)
    sdof_indices = findall(sdofs)

    # Initialize the aux. vectors and matrices based on the storage option:
    if options.storage isa DenseStorage
        K     = zeros(T, num_dofs, num_dofs)
        F_ref = zeros(T, num_dofs)
        F_ext = zeros(T, num_dofs)
        F_int = zeros(T, num_dofs)
        F_res = zeros(T, num_dofs)
        U     = zeros(T, num_dofs)
        δU_p  = zeros(T, num_dofs)
        δU_r  = zeros(T, num_dofs)
        K_red = zeros(T, num_fdofs, num_fdofs)
        U_red = zeros(T, num_fdofs)
        factorization = LUFactorization()
    end

    # Initialize the aux. vectors and matrices based on the storage option:
    if options.storage isa SparseStorage
        K     = spzeros(T, num_dofs, num_dofs)
        F_ref = spzeros(T, num_dofs)
        F_ext = spzeros(T, num_dofs)
        F_int = spzeros(T, num_dofs)
        F_res = spzeros(T, num_dofs)
        U     = spzeros(T, num_dofs)
        δU_p  = spzeros(T, num_dofs)
        δU_r  = spzeros(T, num_dofs)
        K_red = spzeros(T, num_fdofs, num_fdofs)
        U_red = spzeros(T, num_fdofs)
        factorization = UMFPACKFactorization()
    end

    # Initialize the reduced force vector:
    # NOTE: It is initialized as a dense vector regardless of the storage option, since it is used as the right-hand side of the linear system, which is typically solved using a dense solver.
    F_red = zeros(T, num_fdofs)

    # Assemble the global stiffness matrix:
    _assemble_K!(K, model)

    # Assemble the global force vector:
    _assemble_F!(F_ref, model)

    # Assemble the reduced stiffness matrix:
    _assemble_K_red!(K_red, K, fdof_indices)

    # Initialize the linear solver cache:
    lscache = init(LinearProblem(K_red, F_red), factorization)

    # Return the analysis cache:
    return AnalysisCache(fdofs, sdofs, fdof_indices, sdof_indices, F_ref, F_ext, F_int, F_res, U, δU_p, δU_r, K, K_red, F_red, U_red, lscache)
end
