"""
    RootedSpanningForest{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}

Concrete type reprensenting a point process defined on the edges of a `graph` characterizing the uniform distribution on the [spanning forests](https://en.wikipedia.org/wiki/Spanning_tree) of `graph` rooted at `roots`.

It can be viewed as a the product distribution of the uniform distribution on the set of neighbors of each vertex conditioned on forming no cycles.

The object has two fields:

- `graph::LG.SimpleGraph{Int64}`
- `roots::Set{Int64}`

**See also**

Section 4.2 of [GuJeLi19](@cite)
"""
struct RootedSpanningForest{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}
    "Graph"
    graph::LG.SimpleGraph{Int64}
    "Roots"
    roots::Set{Int64}
end

"""
    RootedSpanningForest(
        graph::LG.SimpleGraph{T},
        roots::Union{Nothing,T,AbstractVector{T},AbstractSet{T}}=nothing
    ) where {T<:Integer}

Construct a [`RootedSpanningForest`](@ref) model on `graph`, rooted at `roots`
"""
function RootedSpanningForest(
    graph::LG.SimpleGraph{T},
    roots::Union{Nothing,T,AbstractVector{T},AbstractSet{T}}=nothing
) where {T<:Integer}
    roots_ = Set(roots === nothing ? rand(1:LG.nv(graph)) : roots)
    return RootedSpanningForest{LG.SimpleDiGraph{T}}(graph, roots_)
end

"""
    generate_sample(
        pp::RootedSpanningForest{T};
        rng=-1
    )

Generate an exact sample from the [`PRS.RootedSpanningForest`](@ref).

Default sampler is [`PRS.generate_sample_prs`](@ref).
"""
generate_sample(pp::RootedSpanningForest; rng=-1) = generate_sample_prs(pp; rng=rng)

"""
    generate_sample_prs(
        pp::RootedSpanningForest{T};
        rng=-1
    )::T where {T<:LG.SimpleDiGraph{Int64}}

Generate a rooted spanning forest of `pp.graph`, uniformly at random among all rooted spanning forests rooted at `pp.roots`, using Partial Rejection Sampling (PRS), see Section 4.2 of [GuJeLi19](@cite)
"""
function generate_sample_prs(
    pp::RootedSpanningForest{T};
    rng=-1
)::T where {T<:LG.SimpleDiGraph{Int64}}
    return _generate_sample_rooted_spanning_forest_prs(pp.graph, pp.roots; rng=rng)
end

"""
Generate a rooted spanning forest of `graph`, uniformly at random among all rooted spanning forests rooted at `roots`, using Partial Rejection Sampling (PRS), see Section 4.2 of [GuJeLi19](@cite)
"""
function _generate_sample_rooted_spanning_forest_prs(
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

"""
    random_neighbor_assignment(
        graph::LG.SimpleGraph{T},
        roots;
        rng=-1
    )::LG.SimpleDiGraph{T} where {T}

Return a oriented subgraph of `graph` where each vertex except the `roots` is connected to a unique neighbor, i.e., each vertex has outdegree equal to one.
"""
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