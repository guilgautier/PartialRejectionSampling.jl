using PartialRejectionSampling

using Plots
using GraphPlot, Colors

using LightGraphs
const LG = LightGraphs

using Cairo, Compose

function plot(
    graph::LG.AbstractGraph,
    dims::Vector{Int}=zeros(Int, 2),
    file="";
    kwargs...
)
    if isempty(file)
        if any(dims .== 0)
            p = GraphPlot.gplot(graph;
                    nodelabel=LG.vertices(graph),
                    kwargs...
                    )
            return p
        else
            @assert prod(dims) == LG.nv(graph)
            pos = collect(Iterators.product(1:dims[2], 1:dims[1]))[:]
            locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

            p = GraphPlot.gplot(graph,
                                locs_x,
                                reverse(locs_y);
                                nodelabel=LG.vertices(graph),
                                kwargs...
                )
            return p
        end
    else
        if any(dims .== 0)
            Compose.draw(
                PDF(file, 16cm, 16cm),
                GraphPlot.gplot(graph;
                                nodelabel=LG.vertices(graph),
                                kwargs...
                )
            )
        else
            @assert prod(dims) == LG.nv(graph)
            pos = collect(Iterators.product(1:dims[2], 1:dims[1]))[:]
            locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

            Compose.draw(
                PDF(file, 16cm, 16cm),
                GraphPlot.gplot(graph,
                                locs_x,
                                reverse(locs_y);
                                nodelabel=LG.vertices(graph),
                                kwargs...
                )
            )
        end
    end
end

function plot(
    pp::Ising,
    state,
    width::Int=0,
    height::Int=0;
    file=""
)
    c_nodes = ifelse.(
            state .== 1,
            Colors.colorant"gray",
            Colors.colorant"white")

    plot(pp.graph, [width, height], file; nodefillc=c_nodes)
end

function plot(
    pp::HardCoreGraph,
    state,
    width::Int=0,
    height::Int=0;
    file=""
)
    c_nodes = [Colors.colorant"turquoise" for _ in 1:LG.nv(pp.graph)]
    c_nodes[state] .= Colors.colorant"red"

    plot(pp.graph, [width, height], file; nodefillc=c_nodes)
end

function plot(
    pp::SinkFreeGraph,
    sample,
    width::Int=0,
    height::Int=0;
    file=""
)
    c_nodes = [LG.outdegree(sample, v) == 0 ? colorant"red" : colorant"turquoise"
                 for v in LG.vertices(sample)]

    plot(pp.graph, [width, height], file; nodefillc=c_nodes)
end

function plot(
    pp::RootedSpanningForest,
    sample,
    width::Int=0,
    height::Int=0;
    file=""
)
    c_nodes, col_edges = color_cycles(sample)
    c_nodes[collect(pp.roots)] .= Colors.colorant"orange"

    plot(pp.graph, [width, height], file; nodefillc=c_nodes)
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

    col_edges = fill(c_edge, LG.ne(graph))
    col_edges[edge_idx] .= c_edge_cycle

    return c_nodes, col_edges
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

    col_edges = fill(c_edge, LG.ne(graph))
    col_edges[edge_idx] .= c_edge_sink

    return c_nodes, col_edges
end
