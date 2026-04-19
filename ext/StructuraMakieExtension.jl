module StructuraMakieExtension
using Structura
import Structura: plotundeformedmodel, plotundeformedmodel!
import Structura: plotdeformedmodel, plotdeformedmodel!
import Structura: plotfibersection, plotfibersection!
import Structura: plotquadrature, plotquadrature!
using Makie

function sharedmodelattributes()
    Makie.@DocumentedAttributes begin
        nodecolor = :red
        nodemarker = @inherit marker :circle
        nodemarkersize = @inherit markersize 12
        nodestrokecolor = @inherit strokecolor :black
        nodestrokewidth = @inherit strokewidth 1
        nodeidvisible = true
        nodeidalignment = (:left, :bottom)
        nodeidcolor = @inherit nodecolor :red

        elementcolor = :black
        elementlinewidth = @inherit linewidth 1
        elementlinestyle = @inherit linestyle :solid
        elementidvisible = true
        elementidalignment = (:left, :top)
        elementidcolor = @inherit elementcolor :black
    end
end

Makie.@recipe PlotUndeformedModel (structura, ) begin
    sharedmodelattributes()...
    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(p::PlotUndeformedModel)
    # Extract the model:
    structura = p.structura[]

    # Extract the nodes and elements:
    nodes    = structura.nodes
    elements = structura.elements

    # Plot the elements:
    for (id, element) in elements
        node_i = element.node_i
        node_j = element.node_j
        x_i, y_i = coordinates(node_i)
        x_j, y_j = coordinates(node_j)
        lines!(p, [x_i, x_j], [y_i, y_j],
            color     = p.elementcolor,
            linewidth = p.elementlinewidth,
            linestyle = p.elementlinestyle)
        x_m = (x_i + x_j) / 2
        y_m = (y_i + y_j) / 2
        if p.elementidvisible[]
            text!(p, "E$(id)", position = (x_m, y_m),
                color = p.elementidcolor,
                align = p.elementidalignment)
        end
    end

    # Plot the nodes:
    for (id, node) in nodes
        x, y = coordinates(node)
        scatter!(p, x, y,
            color       = p.nodecolor,
            marker      = p.nodemarker,
            markersize  = p.nodemarkersize,
            strokecolor = p.nodestrokecolor,
            strokewidth = p.nodestrokewidth)
        if p.nodeidvisible[]
            text!(p, "N$(id)", position = (x, y),
                color = p.nodeidcolor,
                align = p.nodeidalignment)
        end
    end

    # Return the plot:
    return p
end

Makie.@recipe PlotDeformedModel (structura, scale) begin
    sharedmodelattributes()...
    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(p::PlotDeformedModel)
    # Extract the model:
    structura = p.structura[]

    # Extract the scale factor:
    scale = p.scale[]

    # Extract the nodes and elements:
    nodes    = structura.nodes
    elements = structura.elements

    # Plot the elements:
    for (id, element) in elements
        node_i = element.node_i
        node_j = element.node_j
        x_i, y_i = coordinates(node_i)
        x_j, y_j = coordinates(node_j)
        u_i = node_i.state.u_c
        u_j = node_j.state.u_c
        x_i_deformed = x_i + scale * u_i[1]
        y_i_deformed = y_i + scale * u_i[2]
        x_j_deformed = x_j + scale * u_j[1]
        y_j_deformed = y_j + scale * u_j[2]
        lines!(p, [x_i_deformed, x_j_deformed], [y_i_deformed, y_j_deformed],
            color     = p.elementcolor,
            linewidth = p.elementlinewidth,
            linestyle = p.elementlinestyle)
        x_m = (x_i_deformed + x_j_deformed) / 2
        y_m = (y_i_deformed + y_j_deformed) / 2
        if p.elementidvisible[]
            text!(p, "E$(id)", position = (x_m, y_m),
                color = p.elementidcolor,
                align = p.elementidalignment)
        end
    end

    # Plot the nodes:
    for (id, node) in nodes
        x, y = coordinates(node)
        u = node.state.u_c
        x_deformed = x + scale * u[1]
        y_deformed = y + scale * u[2]
        scatter!(p, x_deformed, y_deformed,
            color       = p.nodecolor,
            marker      = p.nodemarker,
            markersize  = p.nodemarkersize,
            strokecolor = p.nodestrokecolor,
            strokewidth = p.nodestrokewidth)
        if p.nodeidvisible[]
            text!(p, "N$(id)", position = (x_deformed, y_deformed),
                color = p.nodeidcolor,
                align = p.nodeidalignment)
        end
    end

    # Return the plot:
    return p
end

Makie.@recipe PlotFiberSection (section, ) begin
    fibercolor       = :red
    fibermarker      = @inherit marker :circle
    fibermarkersize  = @inherit markersize 12
    fiberstrokecolor = @inherit strokecolor :black
    fiberstrokewidth = @inherit strokewidth 1
    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(p::PlotFiberSection)
    section = p.section[]

    zs = section.zs
    ys = section.ys

    scatter!(p, zs, ys,
        color       = p.fibercolor,
        marker      = p.fibermarker,
        markersize  = p.fibermarkersize,
        strokecolor = p.fiberstrokecolor,
        strokewidth = p.fiberstrokewidth)

    return p
end

Makie.@recipe PlotQuadrature (quadrature, ) begin
    quadraturecolor       = :black
    quadraturemarker      = @inherit marker :circle
    quadraturemarkercolor = :red
    quadraturemarkersize  = @inherit markersize 12
    quadraturestrokecolor = @inherit strokecolor :black
    quadraturestrokewidth = @inherit strokewidth 1
    Makie.mixin_generic_plot_attributes()...
end

function Makie.plot!(p::PlotQuadrature)
    quadrature = p.quadrature[]

    ξ = Vector(quadrature.ξ)
    w = Vector(quadrature.w)

    scatterlines!(p, ξ, w,
        color       = p.quadraturecolor,
        marker      = p.quadraturemarker,
        markercolor = p.quadraturemarkercolor,
        markersize  = p.quadraturemarkersize,
        strokecolor = p.quadraturestrokecolor,
        strokewidth = p.quadraturestrokewidth)

    return p
end
end
