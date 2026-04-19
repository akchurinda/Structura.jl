function _assemble_K(model::Model)
    # Determine the appropriate type:
    T = Float64
    for element in values(model.elements)
        T = promote_type(T, eltype(_get_k_gbl(element)))
    end

    # Initialize the global stiffness matrix:
    num_dofs = 3 * length(model.nodes)
    K = spzeros(T, num_dofs, num_dofs)

    # Assemble the global stiffness matrix:
    _assemble_K!(K, model)

    # Return the global stiffness matrix:
    return K
end

function _assemble_F_ext(model::Model)
    # Determine the appropriate type:
    T = Float64
    for nload in values(model.nloads)
        T = promote_type(T, eltype(nload.f))
    end
    for eload in values(model.eloads)
        T = promote_type(T, eltype(eload.q))
        T = promote_type(T, typeof(_get_L_0(model.elements[eload.eid].transformation)))
    end

    # Initialize the global external force vector:
    num_dofs = 3 * length(model.nodes)
    F = spzeros(T, num_dofs)
    
    # Assemble the global external force vector:
    _assemble_F!(F, model)

    # Return the global external force vector:
    return F
end

function _assemble_F_int(model::Model)
    # Determine the appropriate type:
    T = Float64
    for element in values(model.elements)
        T = promote_type(T, eltype(_get_q_gbl(element)))
    end

    # Initialize the global internal force vector:
    num_dofs = 3 * length(model.nodes)
    F_int = zeros(T, num_dofs)

    # Assemble the global internal force vector:
    _assemble_Q!(F_int, model)
    
    # Return the global internal force vector:
    return F_int
end

function _assemble_K!(K::AbstractMatrix{<:Real}, model::Model)
    # Zero out the matrix:
    fill!(K, zero(eltype(K)))

    # Assemble the global stiffness matrix:
    for element in values(model.elements)
        k = _get_k_gbl(element)

        dofnumbers_i = _get_dofnumbers(element.node_i)
        dofnumbers_j = _get_dofnumbers(element.node_j)

        @inbounds for a in 1:3
            row_i = dofnumbers_i[a]
            row_j = dofnumbers_j[a]

            for b in 1:3
                col_i = dofnumbers_i[b]
                col_j = dofnumbers_j[b]

                K[row_i, col_i] += k[a, b]
                K[row_i, col_j] += k[a, b + 3]
                K[row_j, col_i] += k[a + 3, b]
                K[row_j, col_j] += k[a + 3, b + 3]
            end
        end
    end

    # Return nothing:
    return nothing
end

function _assemble_F!(F::AbstractVector{<:Real}, model::Model)
    # Zero out the vector:
    fill!(F, zero(eltype(F)))

    # Assemble nodal loads:
    for nload in values(model.nloads)
        # Extract the node ID:
        nid = nload.nid

        # Get the DOF numbers for the node:
        dofnumbers = _get_dofnumbers(model.nodes[nid])
        
        # Assemble the nodal load into the global force vector:
        @inbounds for i in 1:3
            F[dofnumbers[i]] += nload.f[i]
        end
    end

    # Assemble element loads:
    for eload in values(model.eloads)
        # Extract the element:
        element = model.elements[eload.eid]

        # Extract the initial length:
        L_0 = _get_L_0(element.transformation)
        
        # Extract the local-to-global transformation matrix:
        T_2_gbl = _get_T_2_gbl(element.transformation)
        
        # Extract the distributed load components:
        q_x, q_y = eload.q
        
        # Compute the fixed-end forces in local coordinates:
        P = q_x * L_0     / 2
        V = q_y * L_0     / 2
        M = q_y * L_0 ^ 2 / 12
        F_lcl = SVector{6}(+P, +V, +M, +P, +V, -M)
        
        # Transform the fixed-end forces to global coordinates:
        F_gbl = transpose(T_2_gbl) * F_lcl
        
        # Get the DOF numbers for the nodes:
        dofnumbers_i = _get_dofnumbers(element.node_i)
        dofnumbers_j = _get_dofnumbers(element.node_j)
        
        # Assemble the fixed-end forces into the global force vector:
        @inbounds for i in 1:3
            F[dofnumbers_i[i]] += F_gbl[i]
            F[dofnumbers_j[i]] += F_gbl[i + 3]
        end
    end

    # Return nothing:
    return nothing
end

function _assemble_Q!(F_int::AbstractVector{<:Real}, model::Model)
    # Zero out the vector:
    fill!(F_int, zero(eltype(F_int)))

    # Assemble the global internal force vector:
    for element in values(model.elements)
        q = _get_q_gbl(element)

        dofnumbers_i = _get_dofnumbers(element.node_i)
        dofnumbers_j = _get_dofnumbers(element.node_j)

        @inbounds for i in 1:3
            F_int[dofnumbers_i[i]] += q[i]
            F_int[dofnumbers_j[i]] += q[i + 3]
        end
    end

    # Return nothing:
    return nothing
end

function _assemble_U!(U::AbstractVector{<:Real}, model::Model)
    # Zero out the vector:
    fill!(U, zero(eltype(U)))

    # Assemble the global displacement vector:
    for node in values(model.nodes)
        dofnumbers = _get_dofnumbers(node)
        @inbounds for i in 1:3
            U[dofnumbers[i]] = node.state.u_c[i]
        end
    end

    # Return nothing:
    return nothing
end

function _assemble_K_red!(K_red::AbstractMatrix{<:Real}, K::AbstractMatrix{<:Real}, fdof_indices::AbstractVector{Int})
    # Extract the free DOF submatrix from the global stiffness matrix
    @inbounds for (j_local, j_global) in enumerate(fdof_indices)
        for (i_local, i_global) in enumerate(fdof_indices)
            K_red[i_local, j_local] = K[i_global, j_global]
        end
    end
    return nothing
end

function _assemble_F_red!(F_red::AbstractVector{<:Real}, F::AbstractVector{<:Real}, fdof_indices::AbstractVector{Int})
    # Extract the free DOF subvector from the global force vector
    @inbounds for (i_local, i_global) in enumerate(fdof_indices)
        F_red[i_local] = F[i_global]
    end
    return nothing
end
