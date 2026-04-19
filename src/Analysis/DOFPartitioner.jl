function partitiondofs(model::Union{ModelBuilder, Model})
    num_dofs = 3 * length(model.nodes)
    fdofs = BitVector(undef, num_dofs)
    sdofs = BitVector(undef, num_dofs)

    fdofs .= 0
    sdofs .= 0

    # Mark only equation numbers that are used after MP reduction as free.
    for node in values(model.nodes)
        dofnumbers = _get_dofnumbers(node)
        @inbounds for i in 1:3
            fdofs[dofnumbers[i]] = 1
        end
    end

    for constraint in model.spconstraints
        node_id = constraint.node_id
        dofnumbers = _get_dofnumbers(model.nodes[node_id])
        for i in 1:3
            if constraint.supports[i]
                fdofs[dofnumbers[i]] = 0
                sdofs[dofnumbers[i]] = 1
            end
        end
    end

    return fdofs, sdofs
end
