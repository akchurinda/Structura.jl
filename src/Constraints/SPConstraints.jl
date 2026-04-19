struct SPConstraint <: AbstractConstraint
    node_id::Int
    supports::SVector{3, Bool}
end
