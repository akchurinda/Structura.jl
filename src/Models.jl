struct Model{
    NT   <: Node,
    ET   <: AbstractElement,
    SPCT <: SPConstraint,
    MPCT <: MPConstraint,
    NLT  <: NLoad,
    ELT  <: ELoad,
    LPFT <: Real}
    nodes        ::OrderedDict{Int, NT}
    elements     ::OrderedDict{Int, ET}
    spconstraints::Vector{SPCT}
    mpconstraints::Vector{MPCT}
    nloads       ::OrderedDict{Int, NLT}
    eloads       ::OrderedDict{Int, ELT}
    λ            ::Base.RefValue{LPFT}
end

@kwdef mutable struct ModelBuilder{
    NT   <: Node,
    ET   <: AbstractElement,
    MT   <: AbstractMaterial,
    ST   <: AbstractSection,
    QT   <: AbstractQuadrature,
    TT   <: Type{<:AbstractTransformation},
    SPCT <: SPConstraint,
    MPCT <: MPConstraint,
    NLT  <: NLoad,
    ELT  <: ELoad,
    LPFT <: Real}
    nodes          ::OrderedDict{Int, NT}  = OrderedDict{Int, NT}()
    elements       ::OrderedDict{Int, ET}  = OrderedDict{Int, ET}()
    materials      ::OrderedDict{Int, MT}  = OrderedDict{Int, MT}()
    sections       ::OrderedDict{Int, ST}  = OrderedDict{Int, ST}()
    quadratures    ::OrderedDict{Int, QT}  = OrderedDict{Int, QT}()
    transformations::OrderedDict{Int, TT}  = OrderedDict{Int, TT}()
    spconstraints  ::Vector{SPCT}          = Vector{SPCT}()
    mpconstraints  ::Vector{MPCT}          = Vector{MPCT}()
    nloads         ::OrderedDict{Int, NLT} = OrderedDict{Int, NLT}()
    eloads         ::OrderedDict{Int, ELT} = OrderedDict{Int, ELT}()
end

function Model(f::Function)
    # Initialize the model builder:
    builder = ModelBuilder{
        Node,
        AbstractElement,
        AbstractMaterial,
        AbstractSection,
        AbstractQuadrature,
        Type{<:AbstractTransformation},
        SPConstraint,
        MPConstraint,
        NLoad,
        ELoad,
        Real}()
    
    # Construct the model:
    f(builder)

    # Infer concrete container element types from the built model:
    # TODO: This is a bit hacky, I need to revisit this design choice later if it causes issues.
    NT   = isempty(builder.nodes   ) ? Node            : typeof(first(values(builder.nodes   )))
    ET   = isempty(builder.elements) ? AbstractElement : typeof(first(values(builder.elements)))
    SPCT = isempty(builder.spconstraints) ? SPConstraint : typeof(first(builder.spconstraints))
    MPCT = isempty(builder.mpconstraints) ? MPConstraint : typeof(first(builder.mpconstraints))
    NLT  = isempty(builder.nloads) ? NLoad : typeof(first(values(builder.nloads)))
    ELT  = isempty(builder.eloads) ? ELoad : typeof(first(values(builder.eloads)))
    TE   = isempty(builder.elements) ? Float64 : mapreduce(_get_T, promote_type, values(builder.elements))
    TNL  = isempty(builder.nloads)   ? TE : mapreduce(_get_T, promote_type, values(builder.nloads))
    TEL  = isempty(builder.eloads)   ? TE : mapreduce(_get_T, promote_type, values(builder.eloads))
    LPFT = promote_type(TE, TNL, TEL)

    # Construct the model:
    model = Model{NT, ET, SPCT, MPCT, NLT, ELT, LPFT}(
        builder.nodes,
        builder.elements,
        builder.spconstraints,
        builder.mpconstraints,
        builder.nloads,
        builder.eloads,
        Ref{LPFT}(0.0))
    
    # Return the model:
    return model
end

function node!(builder::ModelBuilder, nid::Int, x::Real, y::Real)
    # Check if the node ID is already taken:
    if haskey(builder.nodes, nid)
        throw(ArgumentError("Node with ID $(nid) already exists."))
    end

    # Create the node and add it to the model:
    node = Node(x, y)
    builder.nodes[nid] = node

    # Return nothing:
    return nothing
end


function spconstraint!(builder::ModelBuilder, nid::Int, u_x_supp::Bool, u_y_supp::Bool, θ_z_supp::Bool)
    # Check if the node ID exists:
    if !haskey(builder.nodes, nid)
        throw(ArgumentError("Node with ID $(nid) does not exist."))
    end

    # Create the constraint and add it to the model:
    constraint = SPConstraint(nid, SVector{3}(u_x_supp, u_y_supp, θ_z_supp))
    push!(builder.spconstraints, constraint)

    # Return nothing:
    return nothing
end

function mpconstraint!(builder::ModelBuilder, nid_i::Int, nid_j::Int, u_x_conn::Bool, u_y_conn::Bool, θ_z_conn::Bool)
    # Check if the node (i) ID exists:
    if !haskey(builder.nodes, nid_i)
        throw(ArgumentError("Node with ID $(nid_i) does not exist."))
    end

    # Check if the node (j) ID exists:
    if !haskey(builder.nodes, nid_j)
        throw(ArgumentError("Node with ID $(nid_j) does not exist."))
    end

    # Create the constraint and add it to the model:
    constraint = MPConstraint(nid_i, nid_j, SVector{3}(u_x_conn, u_y_conn, θ_z_conn))
    push!(builder.mpconstraints, constraint)

    # Return nothing:
    return nothing
end

function material!(builder::ModelBuilder, id::Int, type::Type{<:AbstractMaterial}, args...)
    # Check if the material ID is already taken:
    if haskey(builder.materials, id)
        throw(ArgumentError("Material with ID $(id) already exists."))
    end

    # Create the material and add it to the model:
    material = type(args...)
    builder.materials[id] = material

    # Return nothing:
    return nothing
end

function section!(builder::ModelBuilder, sid::Int, type::Type{<:AbstractSection}, args...)
    # Check if the section ID is already taken:
    if haskey(builder.sections, sid)
        throw(ArgumentError("Section with ID $(sid) already exists."))
    end

    # Create the section and add it to the model:
    section = type(args...)
    builder.sections[sid] = section

    # Return nothing:
    return nothing
end

function section!(f::Function, builder::ModelBuilder, sid::Int, ::Type{FiberSection})
    # Check if the section ID is already taken:
    if haskey(builder.sections, sid)
        throw(ArgumentError("Section with ID $(sid) already exists."))
    end

    # Initialize the fiber section builder:
    sectionbuilder = FiberSectionBuilder()

    # Construct the fiber section:
    f(sectionbuilder)

    # Check if the section has any fibers:
    if isempty(sectionbuilder.material_ids)
        throw(ArgumentError("FiberSection requires at least one fiber."))
    end

    # Extract the materials for the fibers:
    _materials = AbstractMaterial[]
    sizehint!(_materials, length(sectionbuilder.material_ids))
    for id in sectionbuilder.material_ids
        if !haskey(builder.materials, id)
            throw(ArgumentError("Material with ID $(id) does not exist."))
        end
        push!(_materials, deepcopy(builder.materials[id]))
    end

    # Determine the promoted scalar type for geometry and material properties:
    T_z = isempty(sectionbuilder.zs) ? Float64 : mapreduce(typeof, promote_type, sectionbuilder.zs)
    T_y = isempty(sectionbuilder.ys) ? Float64 : mapreduce(typeof, promote_type, sectionbuilder.ys)
    T_a = isempty(sectionbuilder.as) ? Float64 : mapreduce(typeof, promote_type, sectionbuilder.as)
    T_m = isempty(_materials)        ? Float64 : mapreduce(_get_T, promote_type, _materials)
    T   = promote_type(T_z, T_y, T_a, T_m)

    MT = typeof(_materials[1])
    if any(!(material isa MT) for material in _materials)
        throw(ArgumentError("FiberSection currently requires all fibers to use the same concrete material type."))
    end

    # Convert the fiber properties to the promoted scalar type:
    zs = convert(Vector{T}, sectionbuilder.zs)
    ys = convert(Vector{T}, sectionbuilder.ys)
    as = convert(Vector{T}, sectionbuilder.as)

    # Compute the centroid of the section and shift the fiber locations:
    a_t = sum(as)
    z_c = sum(z * a for (z, a) in zip(zs, as)) / a_t
    y_c = sum(y * a for (y, a) in zip(ys, as)) / a_t
    zs .-= z_c
    ys .-= y_c

    materials = Vector{MT}(undef, length(_materials))
    for (index, material) in enumerate(_materials)
        materials[index] = material
    end

    # Create the section and add it to the model:
    section = FiberSection{T, MT}(zs, ys, as, materials, FiberSectionState{T}())
    builder.sections[sid] = section

    # Return nothing:
    return nothing
end

function fiber!(builder::FiberSectionBuilder, z::Real, y::Real, a::Real, id::Int)
    # Add the fiber to the section:
    push!(builder.zs, z)
    push!(builder.ys, y)
    push!(builder.as, a)
    push!(builder.material_ids, id)

    # Return nothing:
    return nothing
end

function quadpatch!(builder::FiberSectionBuilder, id::Int,
    z_i::Real, y_i::Real, 
    z_j::Real, y_j::Real,
    z_k::Real, y_k::Real, 
    z_l::Real, y_l::Real,
    num_fibers_ij::Int, 
    num_fibers_jk::Int)
    # Compute the discretization along each dimension:
    dξ = 2 / num_fibers_ij
    dη = 2 / num_fibers_jk

    # Loop over the fibers and add them to the section:
    for i in 1:num_fibers_ij
        for j in 1:num_fibers_jk
            # Compute the location of the fiber in the natural coordinate system:
            ξ_0 = -1 + (i - 1) * dξ
            η_0 = -1 + (j - 1) * dη
            ξ_1 = ξ_0 + dξ
            η_1 = η_0 + dη
            ξ_c = ξ_0 + dξ / 2
            η_c = η_0 + dη / 2
            
            # Interpolate the coordinates of the fiber in the physical coordinate system:
            z_1, y_1 = _interpolate(ξ_0, η_0, z_i, y_i, z_j, y_j, z_k, y_k, z_l, y_l)
            z_2, y_2 = _interpolate(ξ_1, η_0, z_i, y_i, z_j, y_j, z_k, y_k, z_l, y_l)
            z_3, y_3 = _interpolate(ξ_1, η_1, z_i, y_i, z_j, y_j, z_k, y_k, z_l, y_l)
            z_4, y_4 = _interpolate(ξ_0, η_1, z_i, y_i, z_j, y_j, z_k, y_k, z_l, y_l)
            z_c, y_c = _interpolate(ξ_c, η_c, z_i, y_i, z_j, y_j, z_k, y_k, z_l, y_l)
            
            # Compute the area of the fiber using the shoelace formula:
            a = abs(z_1 * y_2 - y_1 * z_2 + z_2 * y_3 - y_2 * z_3 + z_3 * y_4 - y_3 * z_4 + z_4 * y_1 - y_4 * z_1) / 2
            
            # Add the fiber to the section:
            push!(builder.zs, z_c)
            push!(builder.ys, y_c)
            push!(builder.as, a)
            push!(builder.material_ids, id)
        end
    end
    
    # Return nothing:
    return nothing
end

function circpatch!(builder::FiberSectionBuilder, id::Int,
    z_0::Real, 
    y_0::Real,
    r::Real,
    num_fibers_r::Int,
    num_fibers_θ::Int)
    # Compute the discretization along each dimension:
    dr =     r / num_fibers_r
    dθ = 2 * π / num_fibers_θ

    # Precompute some constants:
    constant = (4 * sin(dθ / 2)) / (3 * dθ)

    # Loop over the fibers and add them to the section:
    for i in 1:num_fibers_r
        # Compute the inner and outer radius of this ring of fibers:
        r_i = (i - 1) * dr  # Inner radius of this ring
        r_o = i * dr        # Outer radius of this ring

        # Compute the centroid radius:
        r_c = constant * ((r_o ^ 3 - r_i ^ 3) / (r_o ^ 2 - r_i ^ 2))

        # Compute the area of the fiber:
        a = (r_o ^ 2 - r_i ^ 2) * dθ / 2

        # Loop over the fibers in this ring:
        for j in 1:num_fibers_θ
            θ_mid = (j - 1 / 2) * dθ

            # Compute the coordinates of the fiber:
            z_c = z_0 + r_c * cos(θ_mid)
            y_c = y_0 + r_c * sin(θ_mid)

            # Add the fiber to the section:
            push!(builder.zs, z_c)
            push!(builder.ys, y_c)
            push!(builder.as, a)
            push!(builder.material_ids, id)
        end
    end

    # Return nothing:
    return nothing
end

function isectionpatch!(builder::FiberSectionBuilder, id::Int,
    d::Real, t_w::Real,
    b_f_t::Real, t_f_t::Real,
    b_f_b::Real, t_f_b::Real,
    num_fibers_hor_w ::Int, num_fibers_ver_w::Int,
    num_fibers_hor_tf::Int, num_fibers_ver_tf::Int,
    num_fibers_hor_bf::Int, num_fibers_ver_bf::Int)
    # Construct the web fibers:
    quadpatch!(builder, id,
        -t_w / 2, -d / 2 + t_f_b, 
        +t_w / 2, -d / 2 + t_f_b, 
        +t_w / 2, +d / 2 - t_f_t, 
        -t_w / 2, +d / 2 - t_f_t, 
        num_fibers_hor_w, num_fibers_ver_w)
    
    # Construct the top flange fibers:
    quadpatch!(builder, id,
        -b_f_t / 2, +d / 2 - t_f_t, 
        +b_f_t / 2, +d / 2 - t_f_t, 
        +b_f_t / 2, +d / 2, 
        -b_f_t / 2, +d / 2, 
        num_fibers_hor_tf, num_fibers_ver_tf)

    # Construct the bottom flange fibers:
    quadpatch!(builder, id,
        -b_f_b / 2, -d / 2, 
        +b_f_b / 2, -d / 2, 
        +b_f_b / 2, -d / 2 + t_f_b, 
        -b_f_b / 2, -d / 2 + t_f_b, 
        num_fibers_hor_bf, num_fibers_ver_bf)

    # Return nothing:
    return nothing
end

function isectionpatch!(builder::FiberSectionBuilder, id::Int,
    label::String,
    num_fibers_hor_w ::Int, num_fibers_ver_w::Int,
    num_fibers_hor_tf::Int, num_fibers_ver_tf::Int,
    num_fibers_hor_bf::Int, num_fibers_ver_bf::Int)
    # Make sure the section name is recognized:
    if !haskey(CATALOG, label)
        throw(ArgumentError("I-section with label $(label) is not recognized."))
    end

    # Extract the dimensions from the catalog:
    d, t_w, b_f, t_f = CATALOG[label]

    # Construct the I-section patch:
    isectionpatch!(builder, id,
        d, t_w, b_f, t_f, b_f, t_f,
        num_fibers_hor_w , num_fibers_ver_w,
        num_fibers_hor_tf, num_fibers_ver_tf,
        num_fibers_hor_bf, num_fibers_ver_bf)

    # Return nothing:
    return nothing
end

function quadrature!(builder::ModelBuilder, qid::Int, type::Type{<:AbstractQuadrature}, n::Int)
    # Check if the quadrature ID is already taken:
    if haskey(builder.quadratures, qid)
        throw(ArgumentError("Quadrature with ID $(qid) already exists."))
    end

    # Create the quadrature and add it to the builder:
    quadrature = type(n)
    builder.quadratures[qid] = quadrature

    # Return nothing:
    return nothing
end

function transformation!(builder::ModelBuilder, tid::Int, type::Type{<:AbstractTransformation})
    # Check if the transformation ID is already taken:
    if haskey(builder.transformations, tid)
        throw(ArgumentError("Transformation with ID $(tid) already exists."))
    end

    # Add the transformation type to the builder:
    builder.transformations[tid] = type

    # Return nothing:
    return nothing
end

function element!(builder::ModelBuilder, eid::Int, type::Type{<:AbstractElement}, nid_i::Int, nid_j::Int, sid::Int, qid::Int, tid::Int)
    # Check if the element ID is already taken:
    if haskey(builder.elements, eid)
        throw(ArgumentError("Element with ID $(eid) already exists."))
    end

    # Check if the node (i) ID exists:
    if !haskey(builder.nodes, nid_i)
        throw(ArgumentError("Node with ID $(nid_i) does not exist."))
    end

    # Check if the node (j) ID exists:
    if !haskey(builder.nodes, nid_j)
        throw(ArgumentError("Node with ID $(nid_j) does not exist."))
    end

    # Check if the section ID exists:
    if !haskey(builder.sections, sid)
        throw(ArgumentError("Section with ID $(sid) does not exist."))
    end

    # Check if the quadrature ID exists:
    if !haskey(builder.quadratures, qid)
        throw(ArgumentError("Quadrature with ID $(qid) does not exist."))
    end

    # Check if the transformation ID exists:
    if !haskey(builder.transformations, tid)
        throw(ArgumentError("Transformation with ID $(tid) does not exist."))
    end

    # Extract the nodes:
    node_i = builder.nodes[nid_i]
    node_j = builder.nodes[nid_j]

    # Extract the section:
    section = builder.sections[sid]

    # Extract the quadrature:
    quadrature = builder.quadratures[qid]

    # Extract and create the transformation:
    transformation = builder.transformations[tid](node_i, node_j)

    # Create the element and add it to the model:
    element = type(node_i, node_j, section, quadrature, transformation)
    builder.elements[eid] = element

    # Return nothing:
    return nothing
end

function nload!(builder::ModelBuilder, nid::Int, f_x::Real, f_y::Real, m_z::Real)
    # Check if the node ID exists:
    if !haskey(builder.nodes, nid)
        throw(ArgumentError("Node with ID $(nid) does not exist."))
    end

    # Create the nodal load and add it to the model:
    load = NLoad(nid, SVector{3}(f_x, f_y, m_z))
    builder.nloads[nid] = load

    # Return nothing:
    return nothing
end

function eload!(builder::ModelBuilder, eid::Int, q_x::Real, q_y::Real)
    # Check if the element ID exists:
    if !haskey(builder.elements, eid)
        throw(ArgumentError("Element with ID $(eid) does not exist."))
    end

    # Create the element load and add it to the model:
    load = ELoad(eid, SVector{2}(q_x, q_y))
    builder.eloads[eid] = load

    # Return nothing:
    return nothing
end

function _interpolate(
    ξ::Real, η::Real,
    z_i::Real, y_i::Real, z_j::Real, y_j::Real,
    z_k::Real, y_k::Real, z_l::Real, y_l::Real)
    # Compute the interpolation functions:
    N_1 = (1 - ξ) * (1 - η) / 4
    N_2 = (1 + ξ) * (1 - η) / 4
    N_3 = (1 + ξ) * (1 + η) / 4
    N_4 = (1 - ξ) * (1 + η) / 4
    
    # Interpolate the coordinates:
    z = N_1 * z_i + N_2 * z_j + N_3 * z_k + N_4 * z_l
    y = N_1 * y_i + N_2 * y_j + N_3 * y_k + N_4 * y_l
    
    # Return the interpolated coordinates:
    return z, y
end

function _get_T(model::Model)
    T = Float64

    for element in values(model.elements)
        T = promote_type(T, _get_T(element))
    end

    for nload in values(model.nloads)
        T = promote_type(T, _get_T(nload))
    end

    for eload in values(model.eloads)
        T = promote_type(T, _get_T(eload))
    end

    return T
end
