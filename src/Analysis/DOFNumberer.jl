struct PlainNumberer <: AbstractDOFNumberer end

@kwdef struct CSGNumberer <: AbstractDOFNumberer solver::MatrixBandwidth.Minimization.AbstractSolver = MatrixBandwidth.Minimization.CapraraSalazarGonzalez() end
@kwdef struct DCMNumberer <: AbstractDOFNumberer solver::MatrixBandwidth.Minimization.AbstractSolver = MatrixBandwidth.Minimization.DelCorsoManzini()        end
@kwdef struct SGSNumberer <: AbstractDOFNumberer solver::MatrixBandwidth.Minimization.AbstractSolver = MatrixBandwidth.Minimization.SaxeGurariSudborough()   end
@kwdef struct BFSNumberer <: AbstractDOFNumberer solver::MatrixBandwidth.Minimization.AbstractSolver = MatrixBandwidth.Minimization.BruteForceSearch()       end
@kwdef struct GPSNumberer <: AbstractDOFNumberer solver::MatrixBandwidth.Minimization.AbstractSolver = MatrixBandwidth.Minimization.GibbsPooleStockmeyer()   end
@kwdef struct RCMNumberer <: AbstractDOFNumberer solver::MatrixBandwidth.Minimization.AbstractSolver = MatrixBandwidth.Minimization.ReverseCuthillMcKee()    end

function _adjacency(model::Model)
	nodes = collect(values(model.nodes))
	n = length(nodes)

	node_to_index = IdDict{Any, Int}()
	rows = zeros(Int, n + 2 * length(model.elements) + 2 * length(model.mpconstraints))
	cols = zeros(Int, n + 2 * length(model.elements) + 2 * length(model.mpconstraints))
	vals = zeros(Int, n + 2 * length(model.elements) + 2 * length(model.mpconstraints))

	for (index, node) in enumerate(nodes)
		node_to_index[node] = index
		push!(rows, index)
		push!(cols, index)
		push!(vals, 1)
	end

	for element in values(model.elements)
        node_i = element.node_i
        node_j = element.node_j
		_edge!(rows, cols, vals, node_to_index, node_i, node_j)
	end

	for constraint in model.mpconstraints
        node_i = model.nodes[constraint.node_i_id]
        node_j = model.nodes[constraint.node_j_id]
		_edge!(rows, cols, vals, node_to_index, node_i, node_j)
	end

	return sparse(rows, cols, vals, n, n), nodes
end

function numberdofs!(model::Model, ::PlainNumberer)
    num_dofs = 3 * length(model.nodes)
	if num_dofs == 0
		return nothing
	end

	ordered_nodes = collect(values(model.nodes))

	# Assign a provisional equation number to each nodal DOF in node insertion order.
	dofcounter = 1
	for node in ordered_nodes
		node.state.dofnumbers = SVector{3}(dofcounter, dofcounter + 1, dofcounter + 2)
		dofcounter += 3
	end

	# Merge equation numbers connected by MP (equal-DOF) constraints.
	parent = collect(1:num_dofs)

	for constraint in model.mpconstraints
		dofnumbers_i = _get_dofnumbers(model.nodes[constraint.node_i_id])
		dofnumbers_j = _get_dofnumbers(model.nodes[constraint.node_j_id])

		@inbounds for i in 1:3
			if constraint.supports[i]
				_unionroots!(parent, dofnumbers_i[i], dofnumbers_j[i])
			end
		end
	end

	# Re-number merged equivalence classes into compact equation IDs in node insertion order.
	representative_to_eqid = Dict{Int, Int}()
	eqcounter = 0
	for node in ordered_nodes
		dofnumbers = _get_dofnumbers(node)

		root_1 = _findroot!(parent, dofnumbers[1])
		root_2 = _findroot!(parent, dofnumbers[2])
		root_3 = _findroot!(parent, dofnumbers[3])

		eq_1 = get!(representative_to_eqid, root_1) do
			eqcounter += 1
			eqcounter
		end
		eq_2 = get!(representative_to_eqid, root_2) do
			eqcounter += 1
			eqcounter
		end
		eq_3 = get!(representative_to_eqid, root_3) do
			eqcounter += 1
			eqcounter
		end

		node.state.dofnumbers = SVector{3}(eq_1, eq_2, eq_3)
	end

	return nothing
end

function numberdofs!(model::Model, dofnumberer::AbstractDOFNumberer)
    num_dofs = 3 * length(model.nodes)
	if num_dofs == 0
		return nothing
	end

	adjacency, nodes = _adjacency(model)
	ordering = minimize_bandwidth(adjacency, dofnumberer.solver).ordering
	ordered_nodes = nodes[ordering]

	# Assign a provisional equation number to each nodal DOF in the bandwidth-reduced order.
	dofcounter = 1
	for node in ordered_nodes
		node.state.dofnumbers = SVector{3}(dofcounter, dofcounter + 1, dofcounter + 2)
		dofcounter += 3
	end

	# Merge equation numbers connected by MP (equal-DOF) constraints.
	num_dofs = 3 * length(model.nodes)
	parent = collect(1:num_dofs)

	for constraint in model.mpconstraints
		dofnumbers_i = _get_dofnumbers(model.nodes[constraint.node_i_id])
		dofnumbers_j = _get_dofnumbers(model.nodes[constraint.node_j_id])

		@inbounds for i in 1:3
			if constraint.supports[i]
				_unionroots!(parent, dofnumbers_i[i], dofnumbers_j[i])
			end
		end
	end

	# Re-number merged equivalence classes into compact equation IDs in the bandwidth-reduced order.
	representative_to_eqid = Dict{Int, Int}()
	eqcounter = 0
	for node in ordered_nodes
		dofnumbers = _get_dofnumbers(node)

		root_1 = _findroot!(parent, dofnumbers[1])
		root_2 = _findroot!(parent, dofnumbers[2])
		root_3 = _findroot!(parent, dofnumbers[3])

		eq_1 = get!(representative_to_eqid, root_1) do
			eqcounter += 1
			eqcounter
		end
		eq_2 = get!(representative_to_eqid, root_2) do
			eqcounter += 1
			eqcounter
		end
		eq_3 = get!(representative_to_eqid, root_3) do
			eqcounter += 1
			eqcounter
		end

		node.state.dofnumbers = SVector{3}(eq_1, eq_2, eq_3)
	end

	return nothing
end

function _findroot!(parent::Vector{Int}, i::Int)
    root = i
    while parent[root] != root
        root = parent[root]
    end

    while parent[i] != i
        nxt = parent[i]
        parent[i] = root
        i = nxt
    end

    return root
end

function _unionroots!(parent::Vector{Int}, i::Int, j::Int)
	root_i = _findroot!(parent, i)
	root_j = _findroot!(parent, j)

    if root_i == root_j
        return nothing
    end

    if root_i < root_j
        parent[root_j] = root_i
    else
        parent[root_i] = root_j
    end

    return nothing
end

function _edge!(rows::Vector{Int}, cols::Vector{Int}, vals::Vector{Int}, node_to_index::IdDict{Any, Int}, node_i, node_j)
    index_i = node_to_index[node_i]
    index_j = node_to_index[node_j]

    if index_i == index_j
        return nothing
    end

    push!(rows, index_i)
    push!(cols, index_j)
    push!(vals, 1)
    push!(rows, index_j)
    push!(cols, index_i)
    push!(vals, 1)

    return nothing
end