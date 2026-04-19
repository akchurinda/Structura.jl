mutable struct NodeState{T<:Real}
    # Trial state variables:
    u_t::SVector{3, T} # Global displacements
    q_t::SVector{3, T} # Global internal forces
    r_t::SVector{3, T} # Global residual forces

    # Committed state variables:
    u_c::SVector{3, T} # Global displacements
    q_c::SVector{3, T} # Global internal forces
    r_c::SVector{3, T} # Global residual forces
    
    # Aux. variables:
    dofnumbers::SVector{3, Int} # DOF numbers

    # Inner constructor:
    function NodeState{T}() where {T<:Real}
        # Initialize the static vectors:
        u = SVector{3, T}(zeros(T, 3))
        v = @SVector zeros(Int, 3)

        # Return the node state:
        return new(u, u, u, u, u, u, v)
    end
end

struct Node{T<:Real}
    x::T
    y::T
    state::NodeState
end

Node(x::Real, y::Real) = begin
    T = promote_type(typeof(x), typeof(y))
    Node(promote(x, y)..., NodeState{T}())
end

function _set_trial_state!(node::Node, u_t::SVector{3, <:Real}, q_t::SVector{3, <:Real}, r_t::SVector{3, <:Real})
    # Extract the state:
    state = node.state

    # Update trial state variables:
    state.u_t = u_t
    state.q_t = q_t
    state.r_t = r_t

    # Return nothing:
    return nothing
end

function _commit_state!(node::Node)
    # Extract the state:
    state = node.state

    # Update committed state variables:
    state.u_c = state.u_t
    state.q_c = state.q_t
    state.r_c = state.r_t

    # Return nothing:
    return nothing
end

function _revert_state!(node::Node)
    # Extract the state:
    state = node.state

    # Restore trial state from committed state:
    state.u_t = state.u_c
    state.q_t = state.q_c
    state.r_t = state.r_c

    # Return nothing:
    return nothing
end

coordinates(node::Node) = (node.x, node.y)
_get_u(node::Node) = node.state.u_t
_get_q(node::Node) = node.state.q_t
_get_r(node::Node) = node.state.r_t
_get_dofnumbers(node::Node) = node.state.dofnumbers
