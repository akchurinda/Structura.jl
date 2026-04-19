struct GaussLobattoQuadrature{N} <: AbstractQuadrature{N}
    ξ::SVector{N, Float64}
    w::SVector{N, Float64}
end

function GaussLobattoQuadrature(::Val{N}) where {N}
    ξ, w = gausslobatto(N)
    return GaussLobattoQuadrature{N}(ξ, w)
end

GaussLobattoQuadrature(N::Int) = GaussLobattoQuadrature(Val(N))
