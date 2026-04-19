mutable struct ElasticPerfectlyPlasticMaterialState{T<:Real} <: AbstractMaterialState
    # Trial state variables:
    ϵ_t  ::T
    σ_t  ::T
    E_t  ::T
    ϵ_p_t::T
    
    # Committed state variables:
    ϵ_c  ::T
    σ_c  ::T
    E_c  ::T
    ϵ_p_c::T

    # Inner constructor:
    function ElasticPerfectlyPlasticMaterialState{T}() where {T<:Real}
        return new(zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T))
    end
end

struct ElasticPerfectlyPlasticMaterial{T<:Real} <: AbstractMaterial{T}
    E  ::T
    σ_y::T
    state::ElasticPerfectlyPlasticMaterialState
end

ElasticPerfectlyPlasticMaterial(E::Real, σ_y::Real) = begin
    T = promote_type(typeof(E), typeof(σ_y))
    ElasticPerfectlyPlasticMaterial(promote(E, σ_y)..., ElasticPerfectlyPlasticMaterialState{T}())
end

_get_T(::ElasticPerfectlyPlasticMaterial{T}) where {T <: Real} = T
_get_ϵ(material::ElasticPerfectlyPlasticMaterial) = material.state.ϵ_t
_get_σ(material::ElasticPerfectlyPlasticMaterial) = material.state.σ_t
_get_E(material::ElasticPerfectlyPlasticMaterial) = material.state.E_t

function _set_trial_state!(material::ElasticPerfectlyPlasticMaterial, ϵ_t::T) where {T <: Real}
    # Extract the material properties:
    E   = material.E
    σ_y = material.σ_y

    # Extract the state:
    state = material.state
    
    # Compute the trial stress based on the last commited plastic strain:
    σ_t = E * (ϵ_t - state.ϵ_p_c)

    # Check if the trial state is elastic or plastic:
    if abs(σ_t) <= σ_y
        state.σ_t   = σ_t
        state.ϵ_t   = ϵ_t
        state.E_t   = E
        state.ϵ_p_t = state.ϵ_p_c
    else
        state.σ_t   = σ_y * sign(σ_t)
        state.ϵ_t   = ϵ_t
        state.E_t   = zero(T)
        state.ϵ_p_t = ϵ_t - state.σ_t / E
    end

    # Return nothing:
    return nothing
end

function _commit_state!(material::ElasticPerfectlyPlasticMaterial)
    # Extract the state:
    state = material.state

    # Commit the trial state:
    state.ϵ_c   = state.ϵ_t
    state.σ_c   = state.σ_t
    state.E_c   = state.E_t
    state.ϵ_p_c = state.ϵ_p_t

    # Return nothing:
    return nothing
end

function _revert_state!(material::ElasticPerfectlyPlasticMaterial)
    # Extract the state:
    state = material.state

    # Restore trial state from committed state:
    state.ϵ_t = state.ϵ_c
    state.σ_t = state.σ_c
    state.E_t = state.E_c
    state.ϵ_p_t = state.ϵ_p_c

    # Return nothing:
    return nothing
end

function set_init_strain!(material::ElasticPerfectlyPlasticMaterial, ϵ_init::T) where {T <: Real}
    # Set the trial state:
    _set_trial_state!(material, ϵ_init)

    # Set the committed state:
    _commit_state!(material)

    # Return nothing:
    return nothing
end

function set_init_stress!(material::ElasticPerfectlyPlasticMaterial, σ_init::T) where {T <: Real}
    # Extract the material properties:
    E   = material.E
    σ_y = material.σ_y

    # Check if the initial stress exceeds the yield stress and cap it if necessary:
    if abs(σ_init) > σ_y
        @warn "The initial stress exceeds the yield stress. Capping the initial stress to the yield stress."
        σ_init = σ_y * sign(σ_init)
    end

    # Compute the corresponding initial strain:
    ϵ_init = σ_init / E

    # Set the initial strain:
    set_init_strain!(material, ϵ_init)

    # Return nothing:
    return nothing
end
