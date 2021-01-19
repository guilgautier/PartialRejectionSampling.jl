struct RootedSpanningForest{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}
    "Graph"
    g::LG.SimpleGraph{Int64}
    "Roots"
    roots::Set{Int64}
end

function RootedSpanningForest(
        graph::LG.SimpleGraph{T},
        roots::Union{Nothing,T,AbstractVector{T},AbstractSet{T}}=nothing
) where {T<:Integer}
    roots_ = Set(roots === nothing ? rand(1:LG.nv(graph)) : roots)
    return RootedSpanningForest{LG.SimpleDiGraph{T}}(graph, roots_)
end

"""
Sample rooted spanning forest using Partial Rejection Sampling of [H. Guo, M. Jerrum](https://arxiv.org/pdf/1611.01647.pdf)
"""
function generate_sample_prs(
        rsf::RootedSpanningForest{T};
        rng=-1
)::T where {T<:LG.SimpleDiGraph{Int64}}
    return _generate_sample_rooted_spanning_forest(rsf.g, rsf.roots; rng=rng)
end

function _generate_sample_rooted_spanning_forest(
        graph::LG.SimpleGraph{T},
        roots;
        rng=-1
)::LG.SimpleDiGraph{T} where {T}
    rng = getRNG(rng)
    g = random_neighbor_assignment(graph, roots; rng=rng)
    while true
        vertices_in_cycles = Set(Iterators.flatten(LG.simplecycles(g)))
        isempty(vertices_in_cycles) && break
        # Resample the successor of vertices involved in cycles
        for v in vertices_in_cycles
            # Remove current edge (v, w)
            w = LG.neighbors(g, v)[1]
            LG.rem_edge!(g, v, w)
            # Resample the successor w of v and add edge (v, w)
            w = rand(rng, LG.neighbors(graph, v))
            LG.add_edge!(g, v, w)
        end
    end
    return g
end

function random_neighbor_assignment(
        graph::LG.SimpleGraph{T},
        roots;
        rng=-1
)::LG.SimpleDiGraph{T} where {T}
    rng = getRNG(rng)
    g = LG.SimpleDiGraph(LG.nv(graph))
    for v in LG.vertices(graph)
        if v ∉ roots
            w = rand(rng, LG.neighbors(graph, v))
            LG.add_edge!(g, v, w)
        end
    end
    return g
end

function color_cycles(
        graph::LG.SimpleDiGraph{T}
) where {T}

    edge_map = edgemap(graph)
    edge_idx = T[]
    nodes = Set{T}()
    for cycle in LG.simplecycles(graph)
        union!(nodes, cycle)
        for (x, y) in zip(cycle, circshift(cycle, -1))
            push!(edge_idx, edge_map[LG.Edge(x, y)])
        end
    end

    col_nodes = [colorant"turquoise" for _ in 1:LG.nv(graph)]
    col_nodes[collect(nodes)] .= colorant"red"

    col_edges = [colorant"lightgray" for _ in 1:LG.ne(graph)]
    col_edges[edge_idx] .= colorant"red"

    return col_nodes, col_edges
end
