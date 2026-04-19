struct MPConstraint <: AbstractConstraint
    node_i_id::Int
    node_j_id::Int
    supports::SVector{3, Bool}
end
