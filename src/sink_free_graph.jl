struct SinkFreeGraph{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}
    "Graph"
    g::LG.SimpleGraph{Int64}
end

function SinkFreeGraph(
        graph::LG.SimpleGraph{T}
) where {T<:Integer}
    return SinkFreeGraph{LG.SimpleDiGraph{T}}(graph)
end

"""
Sample uniform sink free orientation of a graph using Partial Rejection Sampling
[H. Guo, M. Jerrum](https://arxiv.org/pdf/1611.01647.pdf)
"""
function generate_sample_prs(
        sfg::SinkFreeGraph;
        rng=-1
)::LG.SimpleDiGraph
    rng = getRNG(rng)
    g = random_edge_orientation(sfg.g; rng=rng)
    while true
        sinks = sink_nodes(g)
        isempty(sinks) && break
        # Resample i.e. flip edges forming a sink with probability 0.5
        for v in sinks
            for w in LG.inneighbors(g, v)
                if rand(rng) < 0.5
                    LG.rem_edge!(g, w, v)
                    LG.add_edge!(g, v, w)
                end
            end
        end
    end
    return g
end

function color_sinks(
        graph::LG.SimpleDiGraph{T}
) where {T}

    nodes = T[]
    edge_idx = T[]
    edge_map = edgemap(graph)
    for v in sink_nodes(graph)
        push!(nodes, v)
        for w in LG.inneighbors(graph, v)
            push!(edge_idx, edge_map[LG.Edge(w, v)])
        end
    end

    col_nodes = [Colors.colorant"turquoise" for _ in 1:LG.nv(graph)]
    col_nodes[nodes] .= Colors.colorant"red"

    col_edges = [Colors.colorant"lightgray" for _ in 1:LG.ne(graph)]
    col_edges[edge_idx] .= Colors.colorant"red"

    return col_nodes, col_edges
end
