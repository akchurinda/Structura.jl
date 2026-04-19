struct GaussLegendreQuadrature{N} <: AbstractQuadrature{N}
    ξ::SVector{N, Float64}
    w::SVector{N, Float64}
end

function GaussLegendreQuadrature(::Val{N}) where {N}
    ξ, w = gausslegendre(N)
    return GaussLegendreQuadrature{N}(ξ, w)
end

GaussLegendreQuadrature(N::Int) = GaussLegendreQuadrature(Val(N))
