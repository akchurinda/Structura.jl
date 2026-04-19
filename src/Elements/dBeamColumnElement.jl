mutable struct dBeamColumnElementState{T<:Real} <: AbstractElementState
    # Trial state variables:
    u_bsc_t::SVector{3, T} # Basic deformations
    q_bsc_t::SVector{3, T} # Basic forces
    k_bsc_t::SMatrix{3, 3, T, 9} # Basic stiffness

    # Committed state variables:
    u_bsc_c::SVector{3, T} # Basic deformations
    q_bsc_c::SVector{3, T} # Basic forces
    k_bsc_c::SMatrix{3, 3, T, 9} # Basic stiffness

    function dBeamColumnElementState{T}() where {T<:Real}
        # Initialize the static vectors:
        v = SVector{3, T}(zeros(T, 3))

        # Initialize the static matrices:
        m = SMatrix{3, 3, T, 9}(zeros(T, 3, 3))

        # Return the element state:
        return new(v, v, m, v, v, m)
    end
end

struct dBeamColumnElement{
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
    state::dBeamColumnElementState
end

function dBeamColumnElement(
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
    state = dBeamColumnElementState{ET}()

    # Create the element:
    element = dBeamColumnElement{ET, NIT, NJT, ST, STS, NIP}(node_i, node_j, sections, quadrature, transformation, state)

    # Set the trial state:
    _set_trial_state!(element, _get_u_gbl(element))

    # Commit the trial state:
    _commit_state!(element)

    # Return the element:
    return element
end

function _set_trial_state!(element::dBeamColumnElement{ET, NIT, NJT, ST, STS, NIP}, u_gbl_t::AbstractVector{T}) where {
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

    # Initialize the trial basic forces and stiffness:
    q_bsc_t = zeros(SVector{3, T})
    k_bsc_t = zeros(SMatrix{3, 3, T})

    # Loop over all sections:
    for i in 1:NIP
        # Extract the integration point and weight:
        ξ = element.quadrature.ξ[i]
        w = element.quadrature.w[i]

        # Map the integration points and weights:
        x = L_0 * (1 + ξ) / 2
        w = L_0 * w / 2

        # Compute the a-matrix at the integration point:
        a = _compute_a(x, L_0)

        # Get the trial state of the section:
        u_s_t = a * u_bsc_t

        # Set the trial state of the section:
        _set_trial_state!(element.sections[i], u_s_t)

        # Get the trial section forces and stiffness:
        q_s_t = _get_q(element.sections[i])
        k_s_t = _get_k(element.sections[i])

        # Compute the contribution of the section to the element's trial basic forces and stiffness:
        q_bsc_t += w * transpose(a) * q_s_t
        k_bsc_t += w * transpose(a) * k_s_t * a
    end

    # Set the trial state of the transformation using the trial basic forces:
    _compute_k_g_lcl_and_try!(element.transformation, q_bsc_t)

    # Set the trial state of the element:
    element.state.u_bsc_t = u_bsc_t
    element.state.q_bsc_t = q_bsc_t
    element.state.k_bsc_t = k_bsc_t

    # Return nothing:
    return nothing
end

function _commit_state!(element::dBeamColumnElement)
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

function _revert_state!(element::dBeamColumnElement)
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

function _get_u_gbl(element::dBeamColumnElement)
    # Extract the global displacements at each node:
    u_gbl_t_i = _get_u(element.node_i)
    u_gbl_t_j = _get_u(element.node_j)

    # Concatenate the global displacements into a single vector:
    u_gbl_t = SVector{6}(u_gbl_t_i..., u_gbl_t_j...)

    # Return the global displacements:
    return u_gbl_t
end

function _get_q_gbl(element::dBeamColumnElement)
    # Compute the global forces from the trial basic forces using the transformation:
    q_bsc_t = element.state.q_bsc_t
    q_gbl_t = _get_q_gbl(element.transformation, q_bsc_t)

    # Return the global forces:
    return q_gbl_t
end

function _get_k_gbl(element::dBeamColumnElement)
    # Compute the global stiffness from the trial basic stiffness using the transformation:
    k_bsc_t = element.state.k_bsc_t
    k_gbl_t = _get_k_t_gbl(element.transformation, k_bsc_t)

    # Return the global stiffness:
    return k_gbl_t
end

function _get_f_gbl(element::dBeamColumnElement)
    # Extract the external forces at each node:
    f_gbl_t_i = get_f(element.node_i)
    f_gbl_t_j = get_f(element.node_j)

    # Concatenate the external forces into a single vector:
    f_gbl_t = SVector{6}(f_gbl_t_i..., f_gbl_t_j...)

    # Return the external forces:
    return f_gbl_t
end

_get_T(::dBeamColumnElement{ET, NIT, NJT, ST, STS, NIP}) where {ET, NIT, NJT, ST, STS, NIP} = ET