struct NLoad{T<:Real}
    nid::Int
    f::SVector{3, T}
end

struct ELoad{T<:Real}
    eid::Int
    q::SVector{2, T}
end

_get_T(::NLoad{T}) where {T <: Real} = T
_get_T(::ELoad{T}) where {T <: Real} = T