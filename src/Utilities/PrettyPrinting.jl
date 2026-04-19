function Base.show(io::IO, model::Model)
    # Extract the nodes and elements for display:
    nodes    = model.nodes
    elements = model.elements

    # Format the output:
    output  = styled""
    output *= styled"{cyan,bold:MODEL}"
    output *= styled"\n├─ $(length(model.nodes)) nodes"
    output *= styled"\n└─ $(length(model.elements)) elements"

    # Print the output:
    print(io, output)
end

function Base.show(io::IO, options::AnalysisOptions)
    # Extract the types of the options for display:
    type_of_storage       = typeof(options.storage)
    type_of_algorithm     = typeof(options.algorithm)
    type_of_controlscheme = typeof(options.controlscheme)
    type_of_testingscheme = typeof(options.testingscheme)
    type_of_dofnumberer   = typeof(options.dofnumberer)

    # Format the output:
    output  = styled""
    output *= styled"{cyan,bold:ANALYSIS OPTIONS}"
    output *= styled"\n├─ Storage:        $(type_of_storage      )"
    output *= styled"\n├─ Algorithm:      $(type_of_algorithm    )"
    output *= styled"\n├─ Control scheme: $(type_of_controlscheme)"
    output *= styled"\n├─ Testing scheme: $(type_of_testingscheme)"
    output *= styled"\n├─ Factorization:  $(type_of_factorization)"
    output *= styled"\n└─ DOF numberer:   $(type_of_dofnumberer  )"

    # Print the output:
    print(io, output)
end

function Base.show(io::IO, ::AnalysisCache)
    # Format the output:
    output  = styled""
    output *= styled"{cyan,bold:ANALYSIS CACHE}"

    # Print the output:
    print(io, output)
end