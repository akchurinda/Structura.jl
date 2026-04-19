function cantilever(P::Real, H::Real, L::Real, E::Real, A::Real, I::Real, element::Type{<:AbstractElement}, transformation::Type{<:AbstractTransformation})
    model = Model() do model
        node!(model,  1, ( 1 - 1) / 10 * L, 0.0)
        node!(model,  2, ( 2 - 1) / 10 * L, 0.0)
        node!(model,  3, ( 3 - 1) / 10 * L, 0.0)
        node!(model,  4, ( 4 - 1) / 10 * L, 0.0)
        node!(model,  5, ( 5 - 1) / 10 * L, 0.0)
        node!(model,  6, ( 6 - 1) / 10 * L, 0.0)
        node!(model,  7, ( 7 - 1) / 10 * L, 0.0)
        node!(model,  8, ( 8 - 1) / 10 * L, 0.0)
        node!(model,  9, ( 9 - 1) / 10 * L, 0.0)
        node!(model, 10, (10 - 1) / 10 * L, 0.0)
        node!(model, 11, (11 - 1) / 10 * L, 0.0)

        spconstraint!(model, 1, true, true, true)

        section!(model, 1, ElasticSection, E, A, I)

        quadrature!(model, 1, GaussLegendreQuadrature, 2)

        transformation!(model, 1, transformation)

        element!(model,  1, element,  1,  2, 1, 1, 1)
        element!(model,  2, element,  2,  3, 1, 1, 1)
        element!(model,  3, element,  3,  4, 1, 1, 1)
        element!(model,  4, element,  4,  5, 1, 1, 1)
        element!(model,  5, element,  5,  6, 1, 1, 1)
        element!(model,  6, element,  6,  7, 1, 1, 1)
        element!(model,  7, element,  7,  8, 1, 1, 1)
        element!(model,  8, element,  8,  9, 1, 1, 1)
        element!(model,  9, element,  9, 10, 1, 1, 1)
        element!(model, 10, element, 10, 11, 1, 1, 1)

        nload!(model, 11, -P, -H, 0.0)
    end

    return model
end

compute_theoretical_1st_order_M_max(x::AbstractVector{<:Real}) = compute_theoretical_1st_order_M_max(x...)
function compute_theoretical_1st_order_M_max(P::Real, H::Real, L::Real, E::Real, A::Real, I::Real)
    M_max = +(H * L)
    return M_max
end

compute_theoretical_1st_order_Δ_max(x::AbstractVector{<:Real}) = compute_theoretical_1st_order_Δ_max(x...)
function compute_theoretical_1st_order_Δ_max(P::Real, H::Real, L::Real, E::Real, A::Real, I::Real)
    Δ_max = -(H * L ^ 3) / (3 * E * I)
    return Δ_max
end

compute_theoretical_2nd_order_M_max(x::AbstractVector{<:Real}) = compute_theoretical_2nd_order_M_max(x...)
function compute_theoretical_2nd_order_M_max(P::Real, H::Real, L::Real, E::Real, A::Real, I::Real)
    α = sqrt((P * L ^ 2) / (E * I))
    M_max = +(H * L) * (tan(α) / α)
    return M_max
end

compute_theoretical_2nd_order_Δ_max(x::AbstractVector{<:Real}) = compute_theoretical_2nd_order_Δ_max(x...)
function compute_theoretical_2nd_order_Δ_max(P::Real, H::Real, L::Real, E::Real, A::Real, I::Real)
    α = sqrt((P * L ^ 2) / (E * I))
    Δ_max = -(H * L ^ 3) / (3 * E * I) * ((3 * (tan(α) - α)) / (α ^ 3))
    return Δ_max
end

@testset "Cantilever Beam: Displacement-Based Beam-Column Element" begin
    # Compute the theoretical results:
    theoretical_1st_order_M_max = compute_theoretical_1st_order_M_max(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0)
    theoretical_1st_order_Δ_max = compute_theoretical_1st_order_Δ_max(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0)
    theoretical_2nd_order_M_max = compute_theoretical_2nd_order_M_max(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0)
    theoretical_2nd_order_Δ_max = compute_theoretical_2nd_order_Δ_max(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0)

    @testset "Cantilever Beam: Linear Elastic Analysis: Load Control" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, dBeamColumnElement, LinearTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, LinearAlgorithm())
            controlscheme!(options, LoadControl(1.0))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        analyze!(model, options, cache, 1)

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_1st_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_1st_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.01)
    end

    @testset "Cantilever Beam: Linear Elastic Analysis: Displacement Control" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, dBeamColumnElement, LinearTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, LinearAlgorithm())
            controlscheme!(options, DisplacementControl(11, 2, theoretical_1st_order_Δ_max))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        analyze!(model, options, cache, 1)

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        @test isapprox(structura_M_max, theoretical_1st_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_1st_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.05)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Load Control: Linear Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, dBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, LinearAlgorithm())
            controlscheme!(options, LoadControl(0.1))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.01)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Load Control: Newton-Raphson Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, dBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, NewtonRaphsonAlgorithm())
            controlscheme!(options, LoadControl(0.1))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.01)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Load Control: Modified Newton-Raphson Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, dBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, ModifiedNewtonRaphsonAlgorithm())
            controlscheme!(options, LoadControl(0.1))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.01)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Displacement Control: Linear Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, dBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, LinearAlgorithm())
            controlscheme!(options, DisplacementControl(11, 2, theoretical_2nd_order_Δ_max / 10))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.05)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Displacement Control: Newton-Raphson Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, dBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, NewtonRaphsonAlgorithm())
            controlscheme!(options, DisplacementControl(11, 2, theoretical_2nd_order_Δ_max / 10))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.05)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Displacement Control: Modified Newton-Raphson Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, dBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, ModifiedNewtonRaphsonAlgorithm())
            controlscheme!(options, DisplacementControl(11, 2, theoretical_2nd_order_Δ_max / 10))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.05)
    end

    @testset "Cantilever Beam: Automatic Differentiation: Linear Elastic Analysis: Load Control" begin
        P = 50.0
        H = 1.0
        L = 180.0
        E = 29000.0
        A = 9.1
        I = 110.0

        function f_M(x::AbstractVector{<:Real})
            P, H, L, E, A, I = x
            Δ_target = compute_theoretical_1st_order_Δ_max(P, H, L, E, A, I)

            model = cantilever(P, H, L, E, A, I, dBeamColumnElement, LinearTransformation)

            options = AnalysisOptions() do options
                storage!(options, DenseStorage())
                algorithm!(options, LinearAlgorithm())
                controlscheme!(options, LoadControl(1.0))
            end

            cache = prepare(model, options)

            analyze!(model, options, cache, 1)

            M_max = get_element_q_gbl(model, 1)[3]

            return M_max
        end

        function f_Δ(x::AbstractVector{<:Real})
            P, H, L, E, A, I = x

            model = cantilever(P, H, L, E, A, I, dBeamColumnElement, LinearTransformation)

            options = AnalysisOptions() do options
                storage!(options, DenseStorage())
                algorithm!(options, LinearAlgorithm())
                controlscheme!(options, LoadControl(1.0))
            end

            cache = prepare(model, options)
            analyze!(model, options, cache, 1)

            Δ_max = get_node_u_gbl(model, 11)[2]

            return Δ_max
        end

        x = [P, H, L, E, A, I]

        theoretical_1st_order_M_max, theoretical_1st_order_∇M_max = value_and_gradient(compute_theoretical_1st_order_M_max, AutoFiniteDiff(), x)
        theoretical_1st_order_Δ_max, theoretical_1st_order_∇Δ_max = value_and_gradient(compute_theoretical_1st_order_Δ_max, AutoFiniteDiff(), x)

        structura_M_max, structura_∇M_max = value_and_gradient(f_M, AutoForwardDiff(), x)
        structura_Δ_max, structura_∇Δ_max = value_and_gradient(f_Δ, AutoForwardDiff(), x)

        @test isapprox(structura_M_max, theoretical_1st_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_1st_order_Δ_max, rtol = 0.01)
        @test isapprox(structura_∇M_max, theoretical_1st_order_∇M_max, rtol = 0.01)
        @test isapprox(structura_∇Δ_max, theoretical_1st_order_∇Δ_max, rtol = 0.01)
    end

    @testset "Cantilever Beam: Automatic Differentiation: Linear Elastic Analysis: Displacement Control" begin
        P = 50.0
        H = 1.0
        L = 180.0
        E = 29000.0
        A = 9.1
        I = 110.0

        function f_M(x::AbstractVector{<:Real})
            P, H, L, E, A, I = x
            Δ_target = compute_theoretical_1st_order_Δ_max(P, H, L, E, A, I)

            model = cantilever(P, H, L, E, A, I, dBeamColumnElement, LinearTransformation)

            options = AnalysisOptions() do options
                storage!(options, DenseStorage())
                algorithm!(options, LinearAlgorithm())
                controlscheme!(options, DisplacementControl(11, 2, Δ_target))
            end

            cache = prepare(model, options)

            analyze!(model, options, cache, 1)

            M_max = get_element_q_gbl(model, 1)[3]

            return M_max
        end

        function f_Δ(x::AbstractVector{<:Real})
            P, H, L, E, A, I = x
            Δ_target = compute_theoretical_1st_order_Δ_max(P, H, L, E, A, I)

            model = cantilever(P, H, L, E, A, I, dBeamColumnElement, LinearTransformation)

            options = AnalysisOptions() do options
                storage!(options, DenseStorage())
                algorithm!(options, LinearAlgorithm())
                controlscheme!(options, DisplacementControl(11, 2, Δ_target))
            end

            cache = prepare(model, options)
            
            analyze!(model, options, cache, 1)

            Δ_max = get_node_u_gbl(model, 11)[2]

            return Δ_max
        end

        x = [P, H, L, E, A, I]

        theoretical_1st_order_M_max, theoretical_1st_order_∇M_max = value_and_gradient(compute_theoretical_1st_order_M_max, AutoFiniteDiff(), x)
        theoretical_1st_order_Δ_max, theoretical_1st_order_∇Δ_max = value_and_gradient(compute_theoretical_1st_order_Δ_max, AutoFiniteDiff(), x)

        structura_M_max, structura_∇M_max = value_and_gradient(f_M, AutoForwardDiff(), x)
        structura_Δ_max, structura_∇Δ_max = value_and_gradient(f_Δ, AutoForwardDiff(), x)

        @test isapprox(structura_M_max, theoretical_1st_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_1st_order_Δ_max, rtol = 0.01)
        @test isapprox(structura_∇M_max, theoretical_1st_order_∇M_max, rtol = 0.01)
        @test isapprox(structura_∇Δ_max, theoretical_1st_order_∇Δ_max, rtol = 0.01)
    end
end

@testset "Cantilever Beam: Force-Based Beam-Column Element" begin
    # Compute the theoretical results:
    theoretical_1st_order_M_max = compute_theoretical_1st_order_M_max(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0)
    theoretical_1st_order_Δ_max = compute_theoretical_1st_order_Δ_max(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0)
    theoretical_2nd_order_M_max = compute_theoretical_2nd_order_M_max(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0)
    theoretical_2nd_order_Δ_max = compute_theoretical_2nd_order_Δ_max(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0)

    @testset "Cantilever Beam: Linear Elastic Analysis: Load Control" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, fBeamColumnElement, LinearTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, LinearAlgorithm())
            controlscheme!(options, LoadControl(1.0))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        analyze!(model, options, cache, 1)

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_1st_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_1st_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.01)
    end

    @testset "Cantilever Beam: Linear Elastic Analysis: Displacement Control" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, fBeamColumnElement, LinearTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, LinearAlgorithm())
            controlscheme!(options, DisplacementControl(11, 2, theoretical_1st_order_Δ_max))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        analyze!(model, options, cache, 1)

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        @test isapprox(structura_M_max, theoretical_1st_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_1st_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.05)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Load Control: Linear Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, fBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, LinearAlgorithm())
            controlscheme!(options, LoadControl(0.1))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.01)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Load Control: Newton-Raphson Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, fBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, NewtonRaphsonAlgorithm())
            controlscheme!(options, LoadControl(0.1))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.01)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Load Control: Modified Newton-Raphson Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, fBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, ModifiedNewtonRaphsonAlgorithm())
            controlscheme!(options, LoadControl(0.1))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.01)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Displacement Control: Linear Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, fBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, LinearAlgorithm())
            controlscheme!(options, DisplacementControl(11, 2, theoretical_2nd_order_Δ_max / 10))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.05)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Displacement Control: Newton-Raphson Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, fBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, NewtonRaphsonAlgorithm())
            controlscheme!(options, DisplacementControl(11, 2, theoretical_2nd_order_Δ_max / 10))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.05)
    end

    @testset "Cantilever Beam: Nonlinear Elastic Analysis: Displacement Control: Modified Newton-Raphson Algorithm" begin
        # Create a cantilever beam model:
        model = cantilever(50.0, 1.0, 180.0, 29000.0, 9.1, 110.0, fBeamColumnElement, CorotationalTransformation)

        # Define the analysis options:
        options = AnalysisOptions() do options
            storage!(options, DenseStorage())
            algorithm!(options, ModifiedNewtonRaphsonAlgorithm())
            controlscheme!(options, DisplacementControl(11, 2, theoretical_2nd_order_Δ_max / 10))
            testingscheme!(options, AbsoluteNormUnbalancedForce(1E-6, 10))
        end

        # Prepare the model:
        cache = prepare(model, options)

        # Run the analysis:
        for _ in 1:10
            analyze!(model, options, cache)
        end

        # Extract the results:
        structura_M_max = get_element_q_gbl(model, 1)[3]
        structura_Δ_max = get_node_u_gbl(model, 11)[2]

        # Compare the results:
        @test isapprox(structura_M_max, theoretical_2nd_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_2nd_order_Δ_max, rtol = 0.01)
        @test isapprox(model.λ[], 1.0, rtol = 0.05)
    end

    @testset "Cantilever Beam: Automatic Differentiation: Linear Elastic Analysis: Load Control" begin
        P = 50.0
        H = 1.0
        L = 180.0
        E = 29000.0
        A = 9.1
        I = 110.0

        function f_M(x::AbstractVector{<:Real})
            P, H, L, E, A, I = x
            Δ_target = compute_theoretical_1st_order_Δ_max(P, H, L, E, A, I)

            model = cantilever(P, H, L, E, A, I, fBeamColumnElement, LinearTransformation)

            options = AnalysisOptions() do options
                storage!(options, DenseStorage())
                algorithm!(options, LinearAlgorithm())
                controlscheme!(options, LoadControl(1.0))
            end

            cache = prepare(model, options)

            analyze!(model, options, cache, 1)

            M_max = get_element_q_gbl(model, 1)[3]

            return M_max
        end

        function f_Δ(x::AbstractVector{<:Real})
            P, H, L, E, A, I = x

            model = cantilever(P, H, L, E, A, I, fBeamColumnElement, LinearTransformation)

            options = AnalysisOptions() do options
                storage!(options, DenseStorage())
                algorithm!(options, LinearAlgorithm())
                controlscheme!(options, LoadControl(1.0))
            end

            cache = prepare(model, options)
            analyze!(model, options, cache, 1)

            Δ_max = get_node_u_gbl(model, 11)[2]

            return Δ_max
        end

        x = [P, H, L, E, A, I]

        theoretical_1st_order_M_max, theoretical_1st_order_∇M_max = value_and_gradient(compute_theoretical_1st_order_M_max, AutoFiniteDiff(), x)
        theoretical_1st_order_Δ_max, theoretical_1st_order_∇Δ_max = value_and_gradient(compute_theoretical_1st_order_Δ_max, AutoFiniteDiff(), x)

        structura_M_max, structura_∇M_max = value_and_gradient(f_M, AutoForwardDiff(), x)
        structura_Δ_max, structura_∇Δ_max = value_and_gradient(f_Δ, AutoForwardDiff(), x)

        @test isapprox(structura_M_max, theoretical_1st_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_1st_order_Δ_max, rtol = 0.01)
        @test isapprox(structura_∇M_max, theoretical_1st_order_∇M_max, rtol = 0.01)
        @test isapprox(structura_∇Δ_max, theoretical_1st_order_∇Δ_max, rtol = 0.01)
    end

    @testset "Cantilever Beam: Automatic Differentiation: Linear Elastic Analysis: Displacement Control" begin
        P = 50.0
        H = 1.0
        L = 180.0
        E = 29000.0
        A = 9.1
        I = 110.0

        function f_M(x::AbstractVector{<:Real})
            P, H, L, E, A, I = x
            Δ_target = compute_theoretical_1st_order_Δ_max(P, H, L, E, A, I)

            model = cantilever(P, H, L, E, A, I, fBeamColumnElement, LinearTransformation)

            options = AnalysisOptions() do options
                storage!(options, DenseStorage())
                algorithm!(options, LinearAlgorithm())
                controlscheme!(options, DisplacementControl(11, 2, Δ_target))
            end

            cache = prepare(model, options)

            analyze!(model, options, cache, 1)

            M_max = get_element_q_gbl(model, 1)[3]

            return M_max
        end

        function f_Δ(x::AbstractVector{<:Real})
            P, H, L, E, A, I = x
            Δ_target = compute_theoretical_1st_order_Δ_max(P, H, L, E, A, I)

            model = cantilever(P, H, L, E, A, I, fBeamColumnElement, LinearTransformation)

            options = AnalysisOptions() do options
                storage!(options, DenseStorage())
                algorithm!(options, LinearAlgorithm())
                controlscheme!(options, DisplacementControl(11, 2, Δ_target))
            end

            cache = prepare(model, options)
            
            analyze!(model, options, cache, 1)

            Δ_max = get_node_u_gbl(model, 11)[2]

            return Δ_max
        end

        x = [P, H, L, E, A, I]

        theoretical_1st_order_M_max, theoretical_1st_order_∇M_max = value_and_gradient(compute_theoretical_1st_order_M_max, AutoFiniteDiff(), x)
        theoretical_1st_order_Δ_max, theoretical_1st_order_∇Δ_max = value_and_gradient(compute_theoretical_1st_order_Δ_max, AutoFiniteDiff(), x)

        structura_M_max, structura_∇M_max = value_and_gradient(f_M, AutoForwardDiff(), x)
        structura_Δ_max, structura_∇Δ_max = value_and_gradient(f_Δ, AutoForwardDiff(), x)

        @test isapprox(structura_M_max, theoretical_1st_order_M_max, rtol = 0.01)
        @test isapprox(structura_Δ_max, theoretical_1st_order_Δ_max, rtol = 0.01)
        @test isapprox(structura_∇M_max, theoretical_1st_order_∇M_max, rtol = 0.01)
        @test isapprox(structura_∇Δ_max, theoretical_1st_order_∇Δ_max, rtol = 0.01)
    end
end
