include("ElasticSection.jl")
include("FiberSection.jl")
include("Common.jl")

const CATALOG = load(joinpath(@__DIR__, "Artifacts", "CATALOG.jld2"), "CATALOG")
