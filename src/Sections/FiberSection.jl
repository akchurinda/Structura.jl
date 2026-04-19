mutable struct FiberSectionState{T<:Real} <: AbstractSectionState
    # Trial state variables:
    u_s_t::SVector{2, T} # Section deformations
    q_s_t::SVector{2, T} # Section forces
    k_s_t::SMatrix{2, 2, T, 4} # Section stiffness

    # Committed state variables:
    u_s_c::SVector{2, T} # Section deformations
    q_s_c::SVector{2, T} # Section forces
    k_s_c::SMatrix{2, 2, T, 4} # Section stiffness

    # Inner constructor:
    function FiberSectionState{T}() where {T<:Real}
        # Initialize the static vectors:
        v = SVector{2, T}(zeros(T, 2))

        # Initialize the static matrices:
        m = SMatrix{2, 2, T, 4}(zeros(T, 2, 2))

        # Return the section state:
        return new(v, v, m, v, v, m)
    end
end

struct FiberSection{T<:Real, MT<:AbstractMaterial{T}} <: AbstractSection{T}
    zs::Vector{T}
    ys::Vector{T}
    as::Vector{T}
    materials::Vector{MT}
    state::FiberSectionState{T}
end

@kwdef struct FiberSectionBuilder
    zs::Vector{<:Real} = Real[]
    ys::Vector{<:Real} = Real[]
    as::Vector{<:Real} = Real[]
    material_ids::Vector{Int} = Int[]
end

_get_T(::FiberSection{T}) where {T <: Real} = T
_get_A(section::FiberSection) = sum(section.as)
_get_I(section::FiberSection) = sum(a * z ^ 2 for (z, a) in zip(section.zs, section.as))
_get_u(section::FiberSection) = section.state.u_s_t
_get_q(section::FiberSection) = section.state.q_s_t
_get_k(section::FiberSection) = section.state.k_s_t

function _set_trial_state!(section::FiberSection, u_s_t::AbstractVector{T}) where {T <: Real}
    # Extract the state:
    state = section.state

    # Extract the trial section deformations:
    ϵ_a_t = u_s_t[1]
    χ_z_t = u_s_t[2]

    # Initialize the section forces and stiffness:
    q_1  = zero(T)
    q_2  = zero(T)
    k_11 = zero(T)
    k_12 = zero(T)
    k_21 = zero(T)
    k_22 = zero(T)

    # Loop over all fibers:
    for (y, a, material) in zip(section.ys, section.as, section.materials)
        # Compute the fiber strain:
        ϵ_f_t = ϵ_a_t - y * χ_z_t
        
        # Set the trial state of the material:
        _set_trial_state!(material, ϵ_f_t)
        
        # Get the trial fiber stiffness and stress:
        E_f_t = _get_E(material)
        σ_f_t = _get_σ(material)
        
        # Compute the stress resultants:
        q_1 += σ_f_t * a
        q_2 -= σ_f_t * a * y

        # Compute the section stiffness:
        k_11 += E_f_t * a
        k_12 -= E_f_t * a * y
        k_21 -= E_f_t * a * y
        k_22 += E_f_t * a * y ^ 2
    end

    q_s_t = SVector{2, T}(q_1, q_2)
    k_s_t = SMatrix{2, 2, T, 4}(k_11, k_12, k_21, k_22)

    # Set the trial state of the section:
    state.u_s_t = u_s_t
    state.q_s_t = q_s_t
    state.k_s_t = k_s_t

    # Return nothing:
    return nothing
end

function _commit_state!(section::FiberSection)
    # Extract the state:
    state = section.state

    # Commit the trial state:
    state.u_s_c = state.u_s_t
    state.q_s_c = state.q_s_t
    state.k_s_c = state.k_s_t

    # Commit the state of each material:
    for material in section.materials
        _commit_state!(material)
    end

    # Return nothing:
    return nothing
end

function _revert_state!(section::FiberSection)
    # Extract the state:
    state = section.state

    # Restore trial state from committed state:
    state.u_s_t = state.u_s_c
    state.q_s_t = state.q_s_c
    state.k_s_t = state.k_s_c

    # Revert all fibers' materials:
    for material in section.materials
        _revert_state!(material)
    end

    # Return nothing:
    return nothing
end
