# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_graph.jl")

dims = [5, 5]
g = LG.grid(dims)

pp = PRS.SinkFreeGraph(g)

using Random
rng = Random.MersenneTwister(123)

c_node_normal = Colors.colorant"turquoise";
c_node_update = Colors.colorant"orange";
c_node_bad = Colors.colorant"red";

c_edge_normal = Colors.colorant"lightgray";
c_edge_update = Colors.colorant"orange";
c_edge_bad = Colors.colorant"red";

edge_map = PRS.edgemap(pp.graph)

path(i) = joinpath("docs/plots/output/sink_free_graph", join([lpad(i, 3, "0"), ".pdf"]))

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
