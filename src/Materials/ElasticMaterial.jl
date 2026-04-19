mutable struct ElasticMaterialState{T<:Real} <: AbstractMaterialState
    # Trial state variables:
    ϵ_t::T
    σ_t::T

    # Committed state variables:
    ϵ_c::T
    σ_c::T

    # Inner constructor:
    function ElasticMaterialState{T}() where {T<:Real}
        return new(zero(T), zero(T), zero(T), zero(T))
    end
end

struct ElasticMaterial{T<:Real} <: AbstractMaterial{T}
    E::T
    state::ElasticMaterialState
end

ElasticMaterial(E::Real) = ElasticMaterial(E, ElasticMaterialState{typeof(E)}())

_get_T(::ElasticMaterial{T}) where {T <: Real} = T
_get_ϵ(material::ElasticMaterial) = material.state.ϵ_t
_get_σ(material::ElasticMaterial) = material.state.σ_t
_get_E(material::ElasticMaterial) = material.E

function _set_trial_state!(material::ElasticMaterial, ϵ_t::T) where {T <: Real}
    # Extract the material properties:
    E = material.E

    # Extract the state:
    state = material.state

    # Set the trial state:
    state.ϵ_t = ϵ_t
    state.σ_t = E * ϵ_t

    # Return nothing:
    return nothing
end

function _commit_state!(material::ElasticMaterial)
    # Extract the state:
    state = material.state

    # Commit the trial state:
    state.ϵ_c = state.ϵ_t
    state.σ_c = state.σ_t

    # Return nothing:
    return nothing
end

function _revert_state!(material::ElasticMaterial)
    # Extract the state:
    state = material.state

    # Restore trial state from committed state:
    state.ϵ_t = state.ϵ_c
    state.σ_t = state.σ_c

    # Return nothing:
    return nothing
end

function set_init_strain!(material::AbstractMaterial, ϵ_init::T) where {T <: Real}
    # Set the trial state:
    _set_trial_state!(material, ϵ_init)

    # Set the committed state:
    _commit_state!(material)

    # Return nothing:
    return nothing
end

function set_init_stress!(material::ElasticMaterial, σ_init::T) where {T <: Real}
    # Extract the material properties:
    E = material.E

    # Compute the corresponding initial strain:
    ϵ_init = σ_init / E

    # Set the initial strain:
    set_init_strain!(material, ϵ_init)

    # Return nothing:
    return nothing
end
