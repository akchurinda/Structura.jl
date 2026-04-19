mutable struct CorotationalTransformationState{T<:Real} <: AbstractTransformationState
    # Trial state variables:
    L_t      ::T
    β_t      ::T
    T_2_lcl_t::SMatrix{3, 6, T, 18}
    k_g_lcl_t::SMatrix{6, 6, T, 36}

    # Committed state variables:
    L_c      ::T
    β_c      ::T
    T_2_lcl_c::SMatrix{3, 6, T, 18}
    k_g_lcl_c::SMatrix{6, 6, T, 36}

    # Constant state variables:
    const L_0::T
    const β_0::T
    const T_2_gbl::SMatrix{6, 6, T, 36}
end

struct CorotationalTransformation <: AbstractTransformation
    state::CorotationalTransformationState
end

function CorotationalTransformation(node_i::Node, node_j::Node)
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
    state = CorotationalTransformationState{T}(
        # Trial state variables:
        L_0,
        β_0,
        T_2_lcl_0,
        k_g_lcl_0,
        # Committed state variables:
        L_0,
        β_0,
        T_2_lcl_0,
        k_g_lcl_0,
        # Constant state variables:
        L_0,
        β_0,
        T_2_gbl_0)

    # Return the transformation:
    return CorotationalTransformation(state)
end

_get_L_0(transformation::CorotationalTransformation) = transformation.state.L_0
_get_β_0(transformation::CorotationalTransformation) = transformation.state.β_0
_get_T_2_lcl(transformation::CorotationalTransformation) = transformation.state.T_2_lcl_t
_get_T_2_gbl(transformation::CorotationalTransformation) = transformation.state.T_2_gbl
_get_k_g_lcl(transformation::CorotationalTransformation) = transformation.state.k_g_lcl_t

function _compute_u_bsc_and_try!(transformation::CorotationalTransformation, u_gbl_t::AbstractVector{<:Real})
    # Extract the transformation state:
    state = transformation.state

    # Compute current length of the element:
    Δx_t = state.L_0 * cos(state.β_0) + u_gbl_t[4] - u_gbl_t[1]
    Δy_t = state.L_0 * sin(state.β_0) + u_gbl_t[5] - u_gbl_t[2]
    L_t  = hypot(Δx_t, Δy_t)

    # Compute the current angle of the element:
    β_t = atan(Δy_t, Δx_t)
    α_t = β_t - state.β_0
    s_α_t = sin(α_t)
    c_α_t = cos(α_t)

    # Compute the types of the state variables:
    T = typeof(state.L_t)
    z = zero(T)
    o = one(T)

    # Compute the basic-to-local transformation matrix:
    T_2_lcl_t = @SMatrix [
        -c_α_t     -s_α_t      z +c_α_t     +s_α_t      z;
        -s_α_t/L_t +c_α_t/L_t +o +s_α_t/L_t -c_α_t/L_t  z;
        -s_α_t/L_t +c_α_t/L_t  z +s_α_t/L_t -c_α_t/L_t +o]

    # Set the trial state:
    state.L_t = L_t
    state.β_t = β_t
    state.T_2_lcl_t = T_2_lcl_t

    # Compute the trial basic deformations from the global displacements:
    u_bsc_t = SVector{3}(L_t - state.L_0, u_gbl_t[3] - α_t, u_gbl_t[6] - α_t)

    # Return the trial basic deformations:
    return u_bsc_t
end

function _compute_k_g_lcl_and_try!(transformation::CorotationalTransformation, q_b_t::AbstractVector{<:Real})
    # Extract the transformation state:
    state = transformation.state

    # Extract the trial state variables:
    L_t = state.L_t
    α_t = state.β_t - state.β_0

    # Extract the trial basic forces:
    N_t, M_1_t, M_2_t = q_b_t[1], q_b_t[2], q_b_t[3]

    # Compute the trigonometric functions of the angle of rotation:
    s_α_t  = sin(α_t)
    c_α_t  = cos(α_t)
    ss_α_t = s_α_t ^ 2
    cc_α_t = c_α_t ^ 2
    sc_α_t = s_α_t * c_α_t

    # Compute the types of the state variables:
    T = typeof(state.L_t)
    z = zero(T)
    o = one(T)

    # Compute the geometric stiffness matrix due to axial forces:
    k_g_a_lcl_t = N_t / L_t * @SMatrix [
        +ss_α_t -sc_α_t z -ss_α_t +sc_α_t z;
        -sc_α_t +cc_α_t z +sc_α_t -cc_α_t z;
              z       z z       z       z z;
        -ss_α_t +sc_α_t z +ss_α_t -sc_α_t z;
        +sc_α_t -cc_α_t z -sc_α_t +cc_α_t z;
              z       z z       z       z z]

    # Compute the geometric stiffness matrix due to bending moments:
    k_g_f_lcl_t = (M_1_t + M_2_t) / L_t ^ 2 * @SMatrix [
             -2*sc_α_t +cc_α_t-ss_α_t z      +2*sc_α_t -cc_α_t+ss_α_t z;
        +cc_α_t-ss_α_t      +2*sc_α_t z -cc_α_t+ss_α_t      -2*sc_α_t z;
                     z              z z              z              z z;
             +2*sc_α_t -cc_α_t+ss_α_t z      -2*sc_α_t +cc_α_t-ss_α_t z;
        -cc_α_t+ss_α_t      -2*sc_α_t z +cc_α_t-ss_α_t      +2*sc_α_t z;
                     z              z z              z              z z]

    # Compute the total geometric stiffness matrix in the local coordinate system:
    k_g_lcl_t = k_g_a_lcl_t + k_g_f_lcl_t

    # Set the trial state:
    state.k_g_lcl_t = k_g_lcl_t

    # Return nothing:
    return nothing
end

function _commit_state!(transformation::CorotationalTransformation)
    # Extract the transformation state:
    state = transformation.state

    # Commit the trial state:
    state.L_c = state.L_t
    state.β_c = state.β_t
    state.T_2_lcl_c = state.T_2_lcl_t
    state.k_g_lcl_c = state.k_g_lcl_t

    # Return nothing:
    return nothing
end

function _revert_state!(transformation::CorotationalTransformation)
    # Extract the transformation state:
    state = transformation.state

    # Restore trial state from committed state:
    state.L_t = state.L_c
    state.β_t = state.β_c
    state.T_2_lcl_t = state.T_2_lcl_c
    state.k_g_lcl_t = state.k_g_lcl_c

    # Return nothing:
    return nothing
end

function _get_q_gbl(transformation::CorotationalTransformation, q_bsc_t::AbstractVector{<:Real})
    # Extract the transformation state:
    state = transformation.state

    # Extract the transformation matrices:
    T_2_lcl = _get_T_2_lcl(transformation)
    T_2_gbl = _get_T_2_gbl(transformation)

    # Compute the global forces from the basic forces:
    q_lcl_t = transpose(T_2_lcl) * q_bsc_t
    q_gbl_t = transpose(T_2_gbl) * q_lcl_t

    # Return the global forces:
    return q_gbl_t
end

function _get_k_t_gbl(transformation::CorotationalTransformation, k_e_bsc_t::AbstractMatrix{<:Real})
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

    # Compute the total global stiffness:
    k_t_gbl_t = k_e_gbl_t + k_g_gbl_t

    # Return the global stiffness:
    return k_t_gbl_t
end
