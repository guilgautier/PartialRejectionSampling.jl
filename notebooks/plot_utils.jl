using PartialRejectionSampling

using Plots
using GraphPlot, Colors

using LightGraphs
const LG = LightGraphs

function plot(
    graph::LG.AbstractGraph,
    dims::Vector{Int}=zeros(Int, 2);
    kwargs...
)
    if any(dims .== 0)
        p = GraphPlot.gplot(graph;
                            nodelabel=LG.vertices(graph), kwargs...)
        return p
    else
        @assert prod(dims) == LG.nv(graph)
        pos = collect(Iterators.product(1:dims[2], 1:dims[1]))[:]
        locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

        p = GraphPlot.gplot(graph, locs_x, reverse(locs_y);
                            nodelabel=LG.vertices(graph), kwargs...)
        return p
    end
end

function plot(
    pp::PRS.Ising,
    state,
    dims=zeros(2);
    kwargs...
)
    c_nodes = ifelse.(
            state .== 1,
            Colors.colorant"gray",
            Colors.colorant"white")

    plot(pp.graph, dims; nodefillc=c_nodes, kwargs...)
end

function plot(
    pp::PRS.HardCoreGraph,
    state,
    dims=zeros(2);
    kwargs...
)
    c_nodes = [Colors.colorant"white" for _ in 1:LG.nv(pp.graph)]
    c_nodes[state] .= Colors.colorant"gray"

    plot(pp.graph, dims; nodefillc=c_nodes, kwargs...)
end

function plot(
    pp::PRS.SinkFreeGraph,
    sample,
    dims=zeros(2);
    kwargs...
)
    c_nodes, c_edges = color_sinks(sample)

    plot(sample, dims; nodefillc=c_nodes, edgestrokec=c_edges, kwargs...)
end

function plot(
    pp::PRS.RootedSpanningForest,
    sample,
    dims=zeros(2);
    kwargs...
)
    c_nodes, c_edges = color_cycles(sample)
    c_nodes[collect(pp.roots)] .= Colors.colorant"lightgreen"

    plot(sample, dims; nodefillc=c_nodes, edgestrokec=c_edges)
end

function color_cycles(
    graph::LG.SimpleDiGraph{T},
    c_node=Colors.colorant"turquoise",
    c_node_cycle=Colors.colorant"red",
    c_edge=Colors.colorant"lightgray",
    c_edge_cycle=Colors.colorant"red"
) where {T}

    edge_map = PRS.edgemap(graph)
    edge_idx = T[]
    nodes = T[]
    for cycle in LG.simplecycles(graph)
        union!(nodes, cycle)
        for (x, y) in zip(cycle, circshift(cycle, -1))
            push!(edge_idx, edge_map[LG.Edge(x, y)])
        end
    end

    c_nodes = fill(c_node, LG.nv(graph))
    c_nodes[nodes] .= c_node_cycle

    c_edges = fill(c_edge, LG.ne(graph))
    c_edges[edge_idx] .= c_edge_cycle

    return c_nodes, c_edges
end

function color_sinks(
    graph::LG.SimpleDiGraph{T},
    c_node=Colors.colorant"turquoise",
    c_node_sink=Colors.colorant"red",
    c_edge=Colors.colorant"lightgray",
    c_edge_sink=Colors.colorant"red"
) where {T}

    nodes = T[]
    edge_idx = T[]
    edge_map = PRS.edgemap(graph)
    for v in PRS.sink_nodes(graph)
        push!(nodes, v)
        for w in LG.inneighbors(graph, v)
            push!(edge_idx, edge_map[LG.Edge(w, v)])
        end
    end

    c_nodes = fill(c_node, LG.nv(graph))
    c_nodes[nodes] .= c_node_sink

    c_edges = fill(c_edge, LG.ne(graph))
    c_edges[edge_idx] .= c_edge_sink

    return c_nodes, c_edges
end

function plot(
    pp::PRS.AbstractSpatialPointProcess,
    points;
    show_center=true,
    radius=0,
    color="white",
    window::Union{Nothing, PRS.AbstractSpatialWindow}=nothing
)
    p = Plots.plot([0], [0],
            label="", legend=false,
            color="white",
            linewidth=0,
            aspect_ratio=:equal,
            grid=false,
            title="")

    for x in (points isa Matrix ? eachcol(points) : points)
        if radius > 0
            plot_disk!(p, x, radius, color)
        end
        if show_center
            Plots.scatter!(p, [x[1]], [x[2]], markersize=2, color=color, grid=false)
        end
    end

    win = isnothing(window) ? PRS.window(pp) : window
    Plots.xlims!(win.c[1], win.c[1] + win.w[1])
    Plots.ylims!(win.c[2], win.c[2] + win.w[win.w isa Real ? 1 : 2])

    return p
end

function plot_disk!(p, center, radius, color=:white)
    θ = LinRange(0, 2π, 15)
    x, y = center[1] .+ radius .* cos.(θ), center[2] .+ radius .* sin.(θ)

    Plots.plot!(p, x, y, seriestype=[:shape], c=color, linewidth=0.5, legend=false, fillapha=0.2)

    return p
end
