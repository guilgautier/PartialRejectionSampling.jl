using PartialRejectionSampling
using Random

using Plots
using GraphPlot, Colors

using LightGraphs
const LG = LightGraphs

using Cairo, Compose

function plot(
    graph::LG.AbstractGraph,
    dims::Vector{Int}=zeros(Int, 2),
    path="";
    kwargs...
)
    if isempty(path)
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
                PDF(path, 16cm, 16cm),
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
                PDF(path, 16cm, 16cm),
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
    pp::PRS.Ising,
    state,
    dims=zeros(2),
    path="";
    kwargs...
)
    c_nodes = ifelse.(
            state .== 1,
            Colors.colorant"gray",
            Colors.colorant"white")

    plot(pp.graph, dims, path; nodefillc=c_nodes, kwargs...)
end

function plot(
    pp::PRS.HardCoreGraph,
    state,
    dims=zeros(2),
    path="";
    kwargs...
)
    c_nodes = [Colors.colorant"white" for _ in 1:LG.nv(pp.graph)]
    c_nodes[state] .= Colors.colorant"grey"

    plot(pp.graph, dims, path; nodefillc=c_nodes, kwargs...)
end

function plot(
    pp::PRS.SinkFreeGraph,
    sample,
    dims=zeros(2),
    path="";
    kwargs...
)
    c_nodes, c_edges = color_sinks(sample)

    plot(sample, dims, path; nodefillc=c_nodes, edgestrokec=c_edges, kwargs...)
end

function plot(
    pp::PRS.RootedSpanningForest,
    sample,
    dims=zeros(2),
    path="";
    kwargs...
)
    c_nodes, c_edges = color_cycles(sample)
    c_nodes[collect(pp.roots)] .= Colors.colorant"orange"

    plot(sample, dims, path; nodefillc=c_nodes, edgestrokec=c_edges)
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

function hard_core_graph_pedagogy(pp, dims, rng)

    path(i) = joinpath(@__DIR__,
                       "output",
                       "hard_core_graph",
                       join([lpad(i, 3, "0"), ".pdf"]))

    c_node_normal = Colors.colorant"turquoise";

    c_node_unccupied = Colors.colorant"white";
    c_node_occupied = Colors.colorant"gray";

    c_node_update = Colors.colorant"orange";
    c_node_bad = Colors.colorant"red";

    i = 0
    c_nodes = fill(c_node_normal, LG.nv(pp.graph))
    plot(pp.graph, dims, path(i); nodefillc=c_nodes)

    proba = pp.β / (one(pp.β) + pp.β)

    adj = LG.adjacency_matrix(pp.graph)
    occupied = Random.randsubseq(rng, LG.vertices(pp.graph), proba)

    while true

        i += 1
        c_nodes = fill(c_node_unccupied, LG.nv(pp.graph))
        c_nodes[occupied] .= c_node_occupied
        plot(pp.graph, dims, path(i); nodefillc=c_nodes)

        # Check if occupied vertices form an independent set
        sub_graph = LG.SimpleGraph(adj[occupied, occupied])
        LG.ne(sub_graph) == 0 && break

        bad = []
        independent, resample = [], Set()
        for cc in LG.connected_components(sub_graph)
            if length(cc) == 1  # Identify current independent vertices
                append!(independent, occupied[cc])
            else  # Construct the resampling set of vertices
                union!(resample, occupied[cc])
                union!(bad, occupied[cc])
                for v in cc
                    union!(resample, LG.neighbors(pp.graph, occupied[v]))
                end
            end
        end

        i += 1
        c_nodes[bad] .= c_node_bad
        plot(pp.graph, dims, path(i); nodefillc=c_nodes)

        Random.randsubseq!(rng, occupied, collect(resample), proba)
        append!(occupied, independent)

        i += 1
        c_nodes[collect(resample)] .= c_node_update;
        plot(pp.graph, dims, path(i); nodefillc=c_nodes)
    end
end

function rooted_spanning_tree_padagogy(pp, dims, rng)

    c_node_roots = Colors.colorant"lightgreen";
    c_node_normal = Colors.colorant"turquoise";
    c_node_update = Colors.colorant"orange";
    c_node_bad = Colors.colorant"red";

    c_edge_normal = Colors.colorant"lightgray";
    c_edge_update = Colors.colorant"orange";
    c_edge_bad = Colors.colorant"red";

    edge_map = PRS.edgemap(pp.graph)

    path(i) = joinpath(@__DIR__,
                       "output",
                       "rooted_spanning_tree",
                       join([lpad(i, 3, "0"), ".pdf"]))

    i = 0
    c_nodes = fill(c_node_normal, LG.nv(pp.graph))
    c_nodes[roots] .= c_node_roots
    plot(pp.graph, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edge_normal)

    i += 1
    g = PRS.random_neighbor_assignment(rng, pp.graph, roots)
    plot(pp.graph, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edge_normal)

    while true
        i += 1
        c_nodes, c_edges = color_cycles(pp.graph, c_node_normal, c_node_bad, c_edge_normal, c_edge_bad)
        c_nodes[roots] .= c_node_roots
        plot(pp.graph, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edges)

        vertices_in_cycles = Set(Iterators.flatten(LG.simplecycles(pp.graph)))
        isempty(vertices_in_cycles) && break

        # Resample the successor of vertices involved in cycles
        edges = []
        c_nodes = fill(c_node_normal, LG.nv(pp.graph))
        c_nodes[roots] .= c_node_roots
        c_edges = fill(c_edge_normal, LG.ne(pp.graph))

        for v in vertices_in_cycles
            c_nodes[v] = c_node_update
            # Remove current edge (v, w)
            w = LG.neighbors(pp.graph, v)[1]
            LG.rem_edge!(pp.graph, v, w)
            # Resample the successor w of v and add edge (v, w)

            w = rand(rng, LG.neighbors(pp.graph, v))
            LG.add_edge!(pp.graph, v, w)

            push!(edges, LG.Edge(v, w))
        end

        i += 1
        edge_map = PRS.edgemap(pp.graph)
        c_edges[[edge_map[e] for e in edges]] .= c_edge_update
        plot(pp.graph, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edges)
    end
end


function sink_free_graph_pedagogy(pp, dims, rng)
    c_node_normal = Colors.colorant"turquoise";
    c_node_update = Colors.colorant"orange";
    c_node_bad = Colors.colorant"red";

    c_edge_normal = Colors.colorant"lightgray";
    c_edge_update = Colors.colorant"orange";
    c_edge_bad = Colors.colorant"red";

    edge_map = PRS.edgemap(pp.graph)

    path(i) = joinpath(@__DIR__,
                       "output",
                       "sink_free_graph",
                       join([lpad(i, 3, "0"), ".pdf"]))

    i = 0
    plot(pp.graph, dims, path(i); nodefillc=c_node_normal, edgestrokec=c_edge_normal)

    i += 1
    g = PRS.random_edge_orientation(rng, pp.graph)
    plot(g, dims, path(i); nodefillc=c_node_normal, edgestrokec=c_edge_normal)

    while true

        i+=1
        c_nodes, c_edges = color_sinks(g, c_node_normal, c_node_bad, c_edge_normal, c_edge_bad)
        plot(g, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edges)

        sinks = PRS.sink_nodes(g)
        isempty(sinks) && break

        # Resample i.e. flip edges forming a sink uniformly at random (Bernoulli(0.5))
        c_nodes = fill(c_node_normal, LG.nv(g))
        c_edges = fill(c_edge_normal, LG.ne(g))

        edges = []
        for v in sinks
            c_nodes[v] = c_node_update
            inneigh_v = copy(LG.inneighbors(g, v))
            for w in inneigh_v
                if rand(rng) < 0.5
                    LG.rem_edge!(g, w, v)
                    LG.add_edge!(g, v, w)
                    push!(edges, LG.Edge(v, w))
                else
                    push!(edges, LG.Edge(w, v))
                end
            end
            edge_map = PRS.edgemap(g)
            c_edges[[edge_map[e] for e in edges]] .= c_edge_update
        end
        i += 1
        plot(g, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edges)
    end
end
