mutable struct ElasticSectionState{T<:Real} <: AbstractSectionState
    # Trial state variables:
    u_s_t::SVector{2, T} # Section deformations
    q_s_t::SVector{2, T} # Section forces

    # Committed state variables:
    u_s_c::SVector{2, T} # Section deformations
    q_s_c::SVector{2, T} # Section forces

    # Constant variables:
    const k_s::SMatrix{2, 2, T, 4} # Section stiffness

    # Inner constructor:
    function ElasticSectionState{T}(E::T, A::T, I::T) where {T<:Real}
        # Initialize the static vectors:
        v = SVector{2, T}(zeros(T, 2))

        # Compute the section stiffness:
        k_s_11 = E * A
        k_s_22 = E * I
        k_s    = @SMatrix T[k_s_11 0; 0 k_s_22]

        # Return the section state:
        return new(v, v, v, v, k_s)
    end
end

struct ElasticSection{T<:Real} <: AbstractSection{T}
    E::T
    A::T
    I::T
    state::ElasticSectionState
end

# Outer constructor for easy initialization
function ElasticSection(E::Real, A::Real, I::Real)
    T = promote_type(typeof(E), typeof(A), typeof(I))
    return ElasticSection(promote(E, A, I)..., ElasticSectionState{T}(promote(E, A, I)...))
end

_get_T(::ElasticSection{T}) where {T <: Real} = T
_get_u(section::ElasticSection) = section.state.u_s_t
_get_q(section::ElasticSection) = section.state.q_s_t
_get_k(section::ElasticSection) = section.state.k_s

function _set_trial_state!(section::ElasticSection, u_s_t::AbstractVector{T}) where {T <: Real}
    # Extract the state:
    state = section.state

    # Compute the section forces:
    q_s_t = state.k_s * u_s_t

    # Set the trial state of the section:
    state.u_s_t = u_s_t
    state.q_s_t = q_s_t

    # Return nothing:
    return nothing
end

function _commit_state!(section::ElasticSection)
    # Extract the state:
    state = section.state

    # Commit the trial state:
    state.u_s_c = state.u_s_t
    state.q_s_c = state.q_s_t

    # Return nothing:
    return nothing
end

function _revert_state!(section::ElasticSection)
    # Extract the state:
    state = section.state

    # Restore trial state from committed state:
    state.u_s_t = state.u_s_c
    state.q_s_t = state.q_s_c

    # Return nothing:
    return nothing
end
