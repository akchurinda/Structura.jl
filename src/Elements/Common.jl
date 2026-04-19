# Define a function that computes the a-matrix (interpolation of deformations):
# NOTE: Used for the displacement-based beam-column element formulation.
function _compute_a(x::Real, L::Real)
    a_11 = 1 / L
    a_12 = 0
    a_13 = 0
    a_21 = 0
    a_22 = 6 * x / L ^ 2 - 4 / L
    a_23 = 6 * x / L ^ 2 - 2 / L
    
    a = SMatrix{2, 3}(a_11, a_21, a_12, a_22, a_13, a_23)

    return a
end

# Define a function that computes the b-matrix (interpolation of forces):
# NOTE: Used for the force-based beam-column element formulation.
function _compute_b(x::Real, L::Real)
    b_11 = 1
    b_12 = 0
    b_13 = 0
    b_21 = 0
    b_22 = x / L - 1
    b_23 = x / L

    b = SMatrix{2, 3}(b_11, b_21, b_12, b_22, b_13, b_23)

    return b
end
