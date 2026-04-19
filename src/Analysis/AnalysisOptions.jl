struct AnalysisOptions
    storage      ::AbstractStorage
    algorithm    ::AbstractAlgorithm    
    controlscheme::AbstractControlScheme
    testingscheme::Union{Nothing, AbstractTestingScheme}
    dofnumberer  ::Union{Nothing, AbstractDOFNumberer  }
end

@kwdef mutable struct AnalysisOptionsBuilder
    storage      ::Union{Nothing, AbstractStorage      } = DenseStorage()
    algorithm    ::Union{Nothing, AbstractAlgorithm    } = nothing
    controlscheme::Union{Nothing, AbstractControlScheme} = nothing
    testingscheme::Union{Nothing, AbstractTestingScheme} = nothing
    dofnumberer  ::Union{Nothing, AbstractDOFNumberer  } = PlainNumberer()
end

function AnalysisOptions(f::Function)
    # Create the builder:
    builder = AnalysisOptionsBuilder()

    # Apply the user function to the builder:
    f(builder)

    # Make sure that the algorithm is defined:
    if isnothing(builder.algorithm)
        throw(ArgumentError("Analysis algorithm is not set in AnalysisOptions."))
    end

    # Make sure that the control scheme is defined:
    if isnothing(builder.controlscheme)
        throw(ArgumentError("Analysis control scheme is not set in AnalysisOptions."))
    end
    
    # Make sure that a testing scheme is provided for the selected algorithm if required:
    if isnothing(builder.testingscheme) && _requires_testing_scheme(builder.algorithm)
        throw(ArgumentError("A testing scheme must be set in AnalysisOptions for nonlinear algorithms."))
    end

    # Create the analysis options:
    options = AnalysisOptions(
        builder.storage,
        builder.algorithm,
        builder.controlscheme,
        builder.testingscheme,
        builder.dofnumberer)

    # Return the analysis options:
    return options
end

storage!(builder::AnalysisOptionsBuilder, storage::Type{<:AbstractStorage}) = storage!(builder, storage())
function storage!(builder::AnalysisOptionsBuilder, storage::AbstractStorage)
    builder.storage = storage
    return nothing
end

algorithm!(builder::AnalysisOptionsBuilder, algorithm::Type{<:AbstractAlgorithm}) = algorithm!(builder, algorithm())
function algorithm!(builder::AnalysisOptionsBuilder, algorithm::AbstractAlgorithm)
    builder.algorithm = algorithm
    return nothing
end

controlscheme!(builder::AnalysisOptionsBuilder, controlscheme::Type{<:AbstractControlScheme}) = controlscheme!(builder, controlscheme())
function controlscheme!(builder::AnalysisOptionsBuilder, controlscheme::AbstractControlScheme)
    builder.controlscheme = controlscheme
    return nothing
end

testingscheme!(builder::AnalysisOptionsBuilder, testingscheme::Type{<:AbstractTestingScheme}) = testingscheme!(builder, testingscheme())
function testingscheme!(builder::AnalysisOptionsBuilder, testingscheme::AbstractTestingScheme)
    builder.testingscheme = testingscheme
    return nothing
end

dofnumberer!(builder::AnalysisOptionsBuilder, dofnumberer::Type{<:AbstractDOFNumberer}) = dofnumberer!(builder, dofnumberer())
function dofnumberer!(builder::AnalysisOptionsBuilder, dofnumberer::AbstractDOFNumberer)
    builder.dofnumberer = dofnumberer
    return nothing
end
