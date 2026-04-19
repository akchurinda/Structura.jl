mutable struct fBeamColumnElementState{T<:Real} <: AbstractElementState
    # Trial state variables:
    u_bsc_t::SVector{3, T} # Basic deformations
    q_bsc_t::SVector{3, T} # Basic forces
    k_bsc_t::SMatrix{3, 3, T, 9} # Basic stiffness

    # Committed state variables:
    u_bsc_c::SVector{3, T} # Basic deformations
    q_bsc_c::SVector{3, T} # Basic forces
    k_bsc_c::SMatrix{3, 3, T, 9} # Basic stiffness

    function fBeamColumnElementState{T}() where {T<:Real}
        # Initialize the static vectors:
        v = SVector{3, T}(zeros(T, 3))

        # Initialize the static matrices:
        m = SMatrix{3, 3, T, 9}(zeros(T, 3, 3))

        # Return the element state:
        return new(v, v, m, v, v, m)
    end
end

struct fBeamColumnElement{
    ET  <: Real,
    NIT <: Real, 
    NJT <: Real, 
    ST  <: Real,
    STS <: AbstractSection{ST}, 
    NIP} <: AbstractElement
    node_i::Node{NIT}
    node_j::Node{NJT}
    sections::SVector{NIP, STS}
    quadrature::AbstractQuadrature{NIP}
    transformation::AbstractTransformation
    state::fBeamColumnElementState
end

function fBeamColumnElement(
    node_i::Node{NIT},
    node_j::Node{NJT},
    section::STS,
    quadrature::AbstractQuadrature{NIP},
    transformation::AbstractTransformation) where {
        NIT <: Real, 
        NJT <: Real,
        ST  <: Real,
        STS <: AbstractSection{ST},
        NIP}
    # Determine the element's scalar type:
    ET = promote_type(NIT, NJT, ST)

    # Create N sections by deep-copying the provided section:
    sections = SVector{NIP, STS}(ntuple(_ -> deepcopy(section), NIP))

    # Initialize the state:
    state = fBeamColumnElementState{ET}()

    # Create the element:
    element = fBeamColumnElement{ET, NIT, NJT, ST, STS, NIP}(node_i, node_j, sections, quadrature, transformation, state)

    # Set the trial state:
    _set_trial_state!(element, _get_u_gbl(element))

    # Commit the trial state:
    _commit_state!(element)

    # Return the element:
    return element
end

const _F_BEAM_COLUMN_ELEMENT_MAXITERATIONS = 25
const _F_BEAM_COLUMN_ELEMENT_TOL = 1E-6
function _set_trial_state!(element::fBeamColumnElement{ET, NIT, NJT, ST, STS, NIP}, u_gbl_t::AbstractVector{T}) where {
    ET  <: Real,
    NIT <: Real,
    NJT <: Real,
    ST  <: Real,
    STS <: AbstractSection{ST},
    NIP,
    T   <: Real}
    # Extract the transformation's properties:
    L_0 = element.transformation.state.L_0
    
    # Compute the trial basic deformations and set the trial state of the transformation:
    u_bsc_t = _compute_u_bsc_and_try!(element.transformation, u_gbl_t)

    # Initialize the trial basic forces and stiffness using the last committed state:
    q_bsc = element.state.q_bsc_c
    k_bsc = element.state.k_bsc_c
    
    # The initial deformation increment to resolve in the loop:
    Δu_bsc = u_bsc_t - element.state.u_bsc_c

    # Element-level state determination loop
    for _ in 1:_F_BEAM_COLUMN_ELEMENT_MAXITERATIONS
        # Update the trial basic forces:
        q_bsc += k_bsc * Δu_bsc
        
        # Initialize the basic element flexibility and residual deformation:
        f_bsc = zeros(SMatrix{3, 3, T})
        u_bsc_res = zeros(SVector{3, T})

        # Loop over all sections:
        for i in 1:NIP
            # Extract the integration point and weight:
            ξ = element.quadrature.ξ[i]
            w = element.quadrature.w[i]

            # Map the integration points and weights:
            x = L_0 * (1 + ξ) / 2
            w = L_0 * w / 2

            # Compute the b-matrix at the integration point:
            b = _compute_b(x, L_0)

            # Target section force based on current element force:
            q_s_target = b * q_bsc
            
            # Retrieve current section trial state
            q_s_prev = _get_q(element.sections[i])
            k_s_prev = _get_k(element.sections[i])
            f_s_prev = inv(k_s_prev) 
            u_s_prev = element.sections[i].state.u_s_t # Extracts current trial deformation
            
            # Section force imbalance:
            Δq_s = q_s_target - q_s_prev
            
            # Estimate section deformation increment and set the trial state:
            Δu_s = f_s_prev * Δq_s
            u_s_trial = u_s_prev + Δu_s
            
            # Set the trial state (Deformation-driven, one pass! No nested loops):
            _set_trial_state!(element.sections[i], u_s_trial)
            
            # Get newly updated section resisting forces and stiffness:
            q_s_resisting = _get_q(element.sections[i])
            k_s = _get_k(element.sections[i])
            f_s = inv(k_s)
            
            # Calculate section force unbalance (Error between target and resisting):
            q_s_unbalanced = q_s_target - q_s_resisting
            
            # Convert force unbalance to a residual section deformation:
            u_s_residual = f_s * q_s_unbalanced
            
            # Integrate element flexibility and residual deformation:
            f_bsc     += w * transpose(b) * f_s * b
            u_bsc_res += w * transpose(b) * (u_s_trial + u_s_residual)
        end

        # Compute the basic stiffness from the basic flexibility:
        k_bsc = inv(f_bsc)
        
        # Error between the integrated section deformations and the target basic deformation:
        Δu_bsc_err = u_bsc_t - u_bsc_res
        
        # Check for convergence:
        if norm(Δu_bsc_err) <= _F_BEAM_COLUMN_ELEMENT_TOL
            break
        end
        
        # If not converged, setup the trial deformation increment for the next iteration:
        Δu_bsc = Δu_bsc_err
    end

    # Set the trial state of the transformation using the trial basic forces:
    _compute_k_g_lcl_and_try!(element.transformation, q_bsc)
    
    # Set the trial state of the element:
    element.state.u_bsc_t = u_bsc_t
    element.state.q_bsc_t = q_bsc
    element.state.k_bsc_t = k_bsc
    
    # Return nothing:
    return nothing
end

function _commit_state!(element::fBeamColumnElement)
    # Commit the element state:
    element.state.u_bsc_c = element.state.u_bsc_t
    element.state.q_bsc_c = element.state.q_bsc_t
    element.state.k_bsc_c = element.state.k_bsc_t

    # Commit the transformation state:
    _commit_state!(element.transformation)

    # Commit the state of each section:
    for section in element.sections
        _commit_state!(section)
    end

    # Return nothing:
    return nothing
end

function _revert_state!(element::fBeamColumnElement)
    # Restore element trial state from committed state:
    element.state.u_bsc_t = element.state.u_bsc_c
    element.state.q_bsc_t = element.state.q_bsc_c
    element.state.k_bsc_t = element.state.k_bsc_c

    # Revert transformation state:
    _revert_state!(element.transformation)

    # Revert all section states:
    for section in element.sections
        _revert_state!(section)
    end

    # Return nothing:
    return nothing
end

function _get_u_gbl(element::fBeamColumnElement)
    # Extract the global displacements at each node:
    u_gbl_t_i = _get_u(element.node_i)
    u_gbl_t_j = _get_u(element.node_j)

    # Concatenate the global displacements into a single vector:
    u_gbl_t = SVector{6}(u_gbl_t_i..., u_gbl_t_j...)

    # Return the global displacements:
    return u_gbl_t
end

function _get_q_gbl(element::fBeamColumnElement)
    # Compute the global forces from the trial basic forces using the transformation:
    q_bsc_t = element.state.q_bsc_t
    q_gbl_t = _get_q_gbl(element.transformation, q_bsc_t)

    # Return the global forces:
    return q_gbl_t
end

function _get_k_gbl(element::fBeamColumnElement)
    # Compute the global stiffness from the trial basic stiffness using the transformation:
    k_bsc_t = element.state.k_bsc_t
    k_gbl_t = _get_k_t_gbl(element.transformation, k_bsc_t)

    # Return the global stiffness:
    return k_gbl_t
end

function _get_f_gbl(element::fBeamColumnElement)
    # Extract the external forces at each node:
    f_gbl_t_i = get_f(element.node_i)
    f_gbl_t_j = get_f(element.node_j)

    # Concatenate the external forces into a single vector:
    f_gbl_t = SVector{6}(f_gbl_t_i..., f_gbl_t_j...)

    # Return the external forces:
    return f_gbl_t
end

_get_T(::fBeamColumnElement{ET, NIT, NJT, ST, STS, NIP}) where {ET, NIT, NJT, ST, STS, NIP} = ET