mutable struct LinearTransformationState{T<:Real} <: AbstractTransformationState
    # Trial state variables:
    # No trial state variables for the linear transformation.

    # Committed state variables:
    # No committed state variables for the linear transformation.
    
    # Constant state variables:
    const L_0    ::T
    const β_0    ::T
    const T_2_lcl::SMatrix{3, 6, T, 18}
    const T_2_gbl::SMatrix{6, 6, T, 36}
    const k_g_lcl::SMatrix{6, 6, T, 36}
end

struct LinearTransformation <: AbstractTransformation
    state::LinearTransformationState
end

function LinearTransformation(node_i::Node, node_j::Node)
    # Compute the initial length of the element:
    Δx_0 = node_j.x - node_i.x
    Δy_0 = node_j.y - node_i.y
    L_0  = hypot(Δx_0, Δy_0)

    # Check that the initial length is larger than zero:
    if iszero(L_0)
        throw(ArgumentError("A transformation requires nodes with nonzero distance."))
    end

    # Compute the initial angle of the element:
    β_0 = atan(Δy_0, Δx_0)
    c_0 = cos(β_0)
    s_0 = sin(β_0)

    # Compute the types of the state variables:
    T = typeof(L_0)
    z = zero(T)
    o = one(T)

    # Compute the local-to-global transformation matrix:
    T_2_gbl_0 = @SMatrix [
        +c_0 +s_0  z    z    z  z;
        -s_0 +c_0  z    z    z  z;
           z    z +o    z    z  z;
           z    z  z +c_0 +s_0  z;
           z    z  z -s_0 +c_0  z;
           z    z  z    z    z +o]

    # Compute the basic-to-local transformation matrix:
    T_2_lcl_0 = @SMatrix [
        -o      z  z +o      z  z;
         z +o/L_0 +o  z -o/L_0  z;
         z +o/L_0  z  z -o/L_0 +o]

    # Compute the geometric stiffness matrix in the local coordinate system:
    k_g_lcl_0 = zero(SMatrix{6, 6, T, 36})

    # Create the transformation state:
    state = LinearTransformationState{T}(L_0, β_0, T_2_lcl_0, T_2_gbl_0, k_g_lcl_0)

    # Return the transformation:
    return LinearTransformation(state)
end

_get_L_0(transformation::LinearTransformation) = transformation.state.L_0
_get_β_0(transformation::LinearTransformation) = transformation.state.β_0
_get_T_2_lcl(transformation::LinearTransformation) = transformation.state.T_2_lcl
_get_T_2_gbl(transformation::LinearTransformation) = transformation.state.T_2_gbl
_get_k_g_lcl(transformation::LinearTransformation) = transformation.state.k_g_lcl

function _compute_u_bsc_and_try!(transformation::LinearTransformation, u_gbl_t::AbstractVector{<:Real})
    # Extract the transformation state:
    state = transformation.state
    
    # Compute the basic deformations from the global displacements:
    u_lcl_t = state.T_2_gbl * u_gbl_t
    u_bsc_t = state.T_2_lcl * u_lcl_t
    
    # Return the trial basic deformations:
    return u_bsc_t
end

function _compute_k_g_lcl_and_try!(transformation::LinearTransformation, q_bsc_t::AbstractVector{<:Real})
    # Return nothing:
    return nothing
end

function _commit_state!(transformation::LinearTransformation)
    # Return nothing:
    return nothing
end

function _revert_state!(transformation::LinearTransformation)
    # Return nothing:
    return nothing
end

function _get_q_gbl(transformation::LinearTransformation, q_bsc_t::AbstractVector{<:Real})
    # Extract the transformation state:
    state = transformation.state

    # Extract the transformation matrices:
    T_2_lcl = _get_T_2_lcl(transformation)
    T_2_gbl = _get_T_2_gbl(transformation)

    # Compute the global forces from the basic forces:
    q_lcl_t = transpose(T_2_lcl) * q_bsc_t
    q_gbl_t = transpose(T_2_gbl) * q_lcl_t
    
    # Return the trial global forces:
    return q_gbl_t
end

function _get_k_t_gbl(transformation::LinearTransformation, k_e_bsc_t::AbstractMatrix{<:Real})
    # Extract the transformation state:
    state = transformation.state

    # Extract the transformation matrices:
    T_2_lcl = _get_T_2_lcl(transformation)
    T_2_gbl = _get_T_2_gbl(transformation)
    
    # Compute the global elastic stiffness from the basic elastic stiffness:
    k_e_lcl_t = transpose(T_2_lcl) * k_e_bsc_t * T_2_lcl
    k_e_gbl_t = transpose(T_2_gbl) * k_e_lcl_t * T_2_gbl

    # Compute the global geometric stiffness from the local geometric stiffness:
    k_g_lcl_t = _get_k_g_lcl(transformation)
    k_g_gbl_t = transpose(T_2_gbl) * k_g_lcl_t * T_2_gbl

    # Compute the global stiffness:
    k_t_gbl_t = k_e_gbl_t + k_g_gbl_t
    
    # Return the trial global stiffness:
    return k_t_gbl_t
end
