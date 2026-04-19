module Structura
# **************************************************
# IMPORTS
# **************************************************
using OrderedCollections
using FastGaussQuadrature
using LinearAlgebra
using StaticArrays
using SparseArrays
using StyledStrings
using LinearSolve
using MatrixBandwidth
using JLD2

# **************************************************
# ABSTRACT COMPONENT STATE TYPES
# **************************************************
abstract type AbstractMaterialState       end
abstract type AbstractSectionState        end
abstract type AbstractTransformationState end
abstract type AbstractElementState        end

# **************************************************
# ABSTRACT COMPONENT TYPES
# **************************************************
abstract type AbstractMaterial{T}    end
abstract type AbstractSection{T}     end
abstract type AbstractQuadrature{N}  end
abstract type AbstractTransformation end
abstract type AbstractElement        end
abstract type AbstractConstraint     end
abstract type AbstractStorage        end
abstract type AbstractAlgorithm      end
abstract type AbstractControlScheme  end
abstract type AbstractTestingScheme  end
abstract type AbstractDOFNumberer    end

# **************************************************
# NODES
# **************************************************
include("Nodes/Nodes.jl")
export NodeState, Node
export coordinates

# **************************************************
# MATERIALS
# **************************************************
include("Materials/Materials.jl")
export AbstractMaterialState, AbstractMaterial
export ElasticMaterialState, ElasticMaterial
export ElasticPerfectlyPlasticMaterialState, ElasticPerfectlyPlasticMaterial

# **************************************************
# SECTIONS
# **************************************************
include("Sections/Sections.jl")
export AbstractSectionState, AbstractSection
export ElasticSectionState, ElasticSection
export FiberSectionState, FiberSection
export CATALOG

# **************************************************
# QUADRATURES
# **************************************************
include("Quadratures/Quadratures.jl")
export AbstractQuadrature
export GaussLegendreQuadrature
export GaussLobattoQuadrature

# **************************************************
# TRANSFORMATIONS
# **************************************************
include("Transformations/Transformations.jl")
export AbstractTransformationState, AbstractTransformation
export LinearTransformationState, LinearTransformation
export CorotationalTransformationState, CorotationalTransformation

# **************************************************
# ELEMENTS
# **************************************************
include("Elements/Elements.jl")
export AbstractElementState, AbstractElement
export dBeamColumnElementState, dBeamColumnElement
export fBeamColumnElementState, fBeamColumnElement

# **************************************************
# CONSTRAINTS
# **************************************************
include("Constraints/Constraints.jl")
export AbstractConstraint
export SPConstraint
export MPConstraint

# **************************************************
# LOADS
# **************************************************
include("Loads/Loads.jl")
export NLoad
export ELoad

# **************************************************
# MODELS
# **************************************************
include("Models.jl")
export Model
export node!
export material!
export section!, fiber!, quadpatch!, circpatch!, isectionpatch!
export quadrature!
export transformation!
export element!
export spconstraint!
export mpconstraint!
export nload!
export eload!

# **************************************************
# ANALYSIS
# **************************************************
include("Analysis/Analysis.jl")
export AbstractStorage
export DenseStorage, SparseStorage
export AbstractAlgorithm
export LinearAlgorithm, NewtonRaphsonAlgorithm, ModifiedNewtonRaphsonAlgorithm
export AbstractControlScheme
export LoadControl, DisplacementControl, WorkControl
export AbstractTestingScheme
export AbsoluteNormUnbalancedForce, AbsoluteNormDisplacementIncrement, AbsoluteNormEnergyIncrement
export AbstractDOFNumberer
export PlainNumberer, CSGNumberer, DCMNumberer, SGSNumberer, BFSNumberer, GPSNumberer, RCMNumberer
export AnalysisOptions
export storage!, algorithm!, controlscheme!, testingscheme!, dofnumberer!
export prepare
export analyze!

# **************************************************
# POST-PROCESSING UTILITIES
# **************************************************
include("Utilities/PostProcessing.jl")
export get_node_u_gbl
export get_node_q_gbl
export get_node_r_gbl
export get_element_u_bsc
export get_element_q_bsc
export get_element_k_bsc
export get_element_u_gbl
export get_element_q_gbl
export get_element_k_gbl
export get_model_λ

# **************************************************
# PLOTTING UTILITIES
# **************************************************
include("Utilities/Plotting.jl")
export plotundeformedmodel
export plotundeformedmodel!
export plotdeformedmodel
export plotdeformedmodel!
export plotfibersection
export plotfibersection!
export plotquadrature
export plotquadrature!

# **************************************************
# PRETTY PRINTING UTILITIES
# **************************************************
include("Utilities/PrettyPrinting.jl")

# **************************************************
# PRECOMPILING UTILITIES
# **************************************************
include("Utilities/Precompiling.jl")
end
