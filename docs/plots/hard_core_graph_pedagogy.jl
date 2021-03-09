# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_graph.jl")
using Random

dims = [5, 5]
g = LG.grid(dims)
β = 0.5

pp = PRS.HardCoreGraph(g, β)

rng = Random.MersenneTwister(123)

path(i) = joinpath("docs/plots/output/hard_core_graph", join([lpad(i, 3, "0"), ".pdf"]))

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
