# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_graph.jl")

dims = [5, 5]
g = LG.grid(dims)
roots = [13]
pp = PRS.RootedSpanningForest(g, roots)

using Random
rng = Random.MersenneTwister(123)

c_node_roots = Colors.colorant"lightgreen";
c_node_normal = Colors.colorant"turquoise";
c_node_update = Colors.colorant"orange";
c_node_bad = Colors.colorant"red";

c_edge_normal = Colors.colorant"lightgray";
c_edge_update = Colors.colorant"orange";
c_edge_bad = Colors.colorant"red";

edge_map = PRS.edgemap(pp.graph)

path(i) = joinpath("docs/plots/output/rooted_spanning_tree", join([lpad(i, 3, "0"), ".pdf"]))

i = 0
c_nodes = fill(c_node_normal, LG.nv(g))
c_nodes[roots] .= c_node_roots
plot(pp.graph, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edge_normal)

i += 1
g = PRS.random_neighbor_assignment(rng, pp.graph, roots)
plot(g, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edge_normal)

while true
    i += 1
    c_nodes, c_edges = color_cycles(g, c_node_normal, c_node_bad, c_edge_normal, c_edge_bad)
    c_nodes[roots] .= c_node_roots
    plot(g, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edges)

    vertices_in_cycles = Set(Iterators.flatten(LG.simplecycles(g)))
    isempty(vertices_in_cycles) && break

    # Resample the successor of vertices involved in cycles
    edges = []
    c_nodes = fill(c_node_normal, LG.nv(g))
    c_nodes[roots] .= c_node_roots
    c_edges = fill(c_edge_normal, LG.ne(g))

    for v in vertices_in_cycles
        c_nodes[v] = c_node_update
        # Remove current edge (v, w)
        w = LG.neighbors(g, v)[1]
        LG.rem_edge!(g, v, w)
        # Resample the successor w of v and add edge (v, w)

        w = rand(rng, LG.neighbors(pp.graph, v))
        LG.add_edge!(g, v, w)

        push!(edges, LG.Edge(v, w))
    end

    i += 1
    edge_map = PRS.edgemap(g)
    c_edges[[edge_map[e] for e in edges]] .= c_edge_update
    plot(g, dims, path(i); nodefillc=c_nodes, edgestrokec=c_edges)
end
