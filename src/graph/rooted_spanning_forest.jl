"""
    RootedSpanningForest{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}

Concrete type reprensenting a point process defined on the edges of a `graph` characterizing the uniform distribution on the [spanning forests](https://en.wikipedia.org/wiki/Spanning_tree) of `graph` rooted at `roots`.

It can be viewed as a the product distribution of the uniform distribution on the set of neighbors of each vertex conditioned on forming no cycles.

The object has two fields:

- `graph::LG.SimpleGraph{Int64}`
- `roots::Vector{Int64}`

**See also**

Section 4.2 of [GuJeLi19](@cite).

# Example

A realization from a ``5\\times 5`` grid graph with `roots=[13]`.

![assets/rooted_spanning_tree.png](assets/rooted_spanning_tree.png)
"""
struct RootedSpanningForest{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}
    "Graph"
    graph::LG.SimpleGraph{Int64}
    "Roots"
    roots::Vector{Int64}
end

function Base.show(io::IO, pp::RootedSpanningForest{T}) where {T}
    print(io, "RootedSpanningForest{$T}\n- graph = $(pp.graph)\n- roots = $(pp.roots)")
end

"""
    RootedSpanningForest(
        graph::LG.SimpleGraph{T},
        roots::Union{Nothing,T,AbstractVector{T},AbstractSet{T}}=nothing
    ) where {T<:Int}

Construct a [`PRS.RootedSpanningForest`](@ref) model on a **connected** `graph`, rooted at `roots`.
If `roots === nothing`, a random vertex is selected uniformly at random among `LG.vertices(g)` and considered as `roots`.

```jldoctest; output = true
using PartialRejectionSampling
using LightGraphs; const LG = LightGraphs

g, roots = LG.grid([5, 5]), [1, 2, 3]

rsf = PRS.RootedSpanningForest(g, roots)

# output

RootedSpanningForest{LightGraphs.SimpleGraphs.SimpleDiGraph{Int64}}
- graph = {25, 40} undirected simple Int64 graph
- roots = [1, 2, 3]
```
"""
function RootedSpanningForest(
    graph::LG.SimpleGraph{T},
    roots::Union{Nothing,T,AbstractVector{T},AbstractSet{T}}=nothing
) where {T<:Int}
    roots_ = unique(roots === nothing ? rand(LG.vertices(graph)) : roots)
    if issubset(roots_, LG.vertices(graph))
        if LG.is_connected(graph)
            return RootedSpanningForest{LG.SimpleDiGraph{T}}(graph, roots_)
        else
            throw(DomainError(graph, "graph should be connected"))
        end
    else
        throw(DomainError(roots, "some roots not contained in vertices(graph)"))
    end
end

"""
    generate_sample(
        [rng::Random.AbstractRNG,]
        pp::RootedSpanningForest
    )

Generate an exact sample from the [`PRS.RootedSpanningForest`](@ref).

Default sampler is [`PRS.generate_sample_prs`](@ref).
"""
function generate_sample(
    rng::Random.AbstractRNG,
    pp::RootedSpanningForest
)
    return generate_sample_prs(rng, pp)
end

# generate_sample(pp::RootedSpanningForest) = generate_sample(Random.default_rng(), pp)

"""
    generate_sample_prs(
        [rng::Random.AbstractRNG,]
        pp::RootedSpanningForest{T}
    )::T where {T<:LG.SimpleDiGraph{Int64}}

Generate a rooted spanning forest of `pp.graph`, uniformly at random among all rooted spanning forests rooted at `pp.roots`, using Partial Rejection Sampling (PRS), see Section 4.2 of [GuJeLi19](@cite).

# Example

An illustration of the procedure on a ``5\\times 5`` grid graph, with `roots=[13]`

- red: variables involved in constraints violation (edges forming a cycle)
- orange: variables resampled (edges originating from a red/orange node)

![assets/rooted_spanning_tree_prs.gif](assets/rooted_spanning_tree_prs.gif)
"""
function generate_sample_prs(
    rng::Random.AbstractRNG,
    pp::RootedSpanningForest{T}
)::T where {T<:LG.SimpleDiGraph{Int64}}
    return _generate_sample_rooted_spanning_forest_prs(rng, pp.graph, pp.roots)
end

# generate_sample_prs(pp::RootedSpanningForest) = generate_sample_prs(Random.default_rng(), pp)

"""
    _generate_sample_rooted_spanning_forest_prs(
        rng::Random.AbstractRNG,
        graph::LG.SimpleGraph{T},
        roots
    )::LG.SimpleDiGraph{T} where {T}

Generate a rooted spanning forest from a **connected** `graph`, uniformly at random among all rooted spanning forests rooted at `roots`, using Wilson's algorithm.
"""
function _generate_sample_rooted_spanning_forest_prs(
    rng::Random.AbstractRNG,
    graph::LG.SimpleGraph{T},
    roots
)::LG.SimpleDiGraph{T} where {T}
    g = random_neighbor_assignment(rng, graph, roots)
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
