@doc raw"""
    sigmoid(x)

Elementwise evaluation of the [sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function)

```math
    \sigma(t) = \frac{1}{1+\exp(t)}
```
"""
sigmoid(x) = @. 1 / (1 + exp(-x))

"""
    allequal(x)::Bool

Inspired by a [Stackoverflow post](https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal/47578613)
"""
@inline function allequal(x)::Bool
    length(x) < 2 && return true
    x1 = x[1]
    @inbounds for i in eachindex(x)
        x[i] == x1 || return false
    end
    return true
end

"""
    normalize_columns!(X::Matrix, p::Real=2)

using `LA.normalize!(x, p)` for each column `x` of `X`.
"""
function normalize_columns!(X::Matrix, p::Real=2)
    for x in eachcol(X)
        LA.normalize!(x, p)
    end
    return X
end

"""
    pairwise_distances(X, Y) = Distances.pairwise(Distances.Euclidean(1e-8), X, Y; dims=2)

Pairwise euclidean distance matrix between columns of `X` and `Y`.
Equivalent to `[norm(x - y) for x in X, y in Y]`.
"""
pairwise_distances(X, Y) = Distances.pairwise(Distances.Euclidean(1e-8), X, Y; dims=2)

"""
    pairwise_distances(X; diag_coeff=Inf) = Distances.pairwise(Distances.Euclidean(1e-8), X; dims=2)

Pairwise euclidean distance matrix between columns of `X`
Equivalent to [`PRS.pairwise_distances`](@ref)`(X, X)` and setting the diagonal elements to `diag_coeff`.
"""
function pairwise_distances(X; diag_coeff=Inf)
    dist = Distances.pairwise(Distances.Euclidean(1e-8), X; dims=2)
    dist[LA.diagind(dist)] .= diag_coeff
    return dist
end

## Graph functions

edgemap(g::LG.AbstractGraph) = Dict(LG.edges(g) .=> 1:LG.ne(g))

sink_nodes(g::LG.SimpleDiGraph) = [v for v in LG.vertices(g) if LG.outdegree(g, v) == 0]

"""
    strong_product(g::G, h::G)::G where {G<:LG.AbstractGraph}

Return the [strong product](https://en.wikipedia.org/wiki/Strong_product_of_graphs) of two graphs

- Section 5 of [BNBLPR19](@cite)
"""
function strong_product(g::G, h::G)::G where {G<:LG.AbstractGraph}
    sp_gh = G(LG.nv(g) * LG.nv(h))
    id(i, j) = (i - 1) * LG.nv(h) + j
    for e1 in LG.edges(g)
        i1, i2 = Tuple(e1)
        for j = 1:LG.nv(h)
            LG.add_edge!(sp_gh, id(i1, j), id(i2, j))
        end
        for e2 in LG.edges(h)
            j1, j2 = Tuple(e2)
            LG.add_edge!(sp_gh, id(i1, j1), id(i2, j2))
            LG.add_edge!(sp_gh, id(i1, j2), id(i2, j1))
        end
    end
    for e in LG.edges(h)
        j1, j2 = Tuple(e)
        for i in LG.vertices(g)
            LG.add_edge!(sp_gh, id(i, j1), id(i, j2))
        end
    end
    return sp_gh
end

"""
    strong_product(g::G)::G where {G<:LG.AbstractGraph}

Return the [strong product](https://en.wikipedia.org/wiki/Strong_product_of_graphs) of a graph with itself

- Section 5 of [BNBLPR19](@cite)
"""
function strong_product(g::G)::G where {G<:LG.AbstractGraph}
    n = LG.nv(g)
    sp_gg = G(n * n)  # nv = n^2, ne = 2n(2n+1)
    id(i, j) = (i - 1) * n + j
    for e1 in LG.edges(g)
        i1, i2 = Tuple(e1)
        for j in LG.vertices(g)
            LG.add_edge!(sp_gg, id(i1, j), id(i2, j))
            LG.add_edge!(sp_gg, id(j, i1), id(j, i2))
        end
        for e2 in LG.edges(g)
            j1, j2 = Tuple(e2)
            LG.add_edge!(sp_gg, id(i1, j1), id(i2, j2))
            LG.add_edge!(sp_gg, id(i1, j2), id(i2, j1))
        end
    end
    return sp_gg
end

king_graph(k::Integer, l::Integer) = strong_product(LG.path_graph(l), LG.path_graph(k))
king_graph(k::Integer) = strong_product(LG.path_graph(k))

"""
    uniform_weighted_graph(
        [rng::Random.AbstractRNG,]
        graph::LG.AbstractGraph
    )::SWG.SimpleWeightedGraph

Return a weighted version of `graph` where each edge is attached an independent uniform random variable.
"""
function uniform_weighted_graph(
    rng::Random.AbstractRNG,
    graph::LG.AbstractGraph
)::SWG.SimpleWeightedGraph
    g = SWG.SimpleWeightedGraph(graph)
    for e in LG.edges(g)
        i, j = Tuple(e)
        @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng, weighttype(g))
    end
    return g
end

uniform_weighted_graph(graph) = uniform_weighted_graph(Random.default_rng(), graph)

"""
    random_edge_orientation(
        rng::Random.AbstractRNG,
        graph::LG.SimpleGraph{T},
        p::Real=0.5
    )::LG.SimpleDiGraph{T} where {T}

Construct an directed version of ``
"""
function random_edge_orientation(
    rng::Random.AbstractRNG,
    graph::LG.SimpleGraph{T},
    p::Real=0.5
)::LG.SimpleDiGraph{T} where {T}
    @assert 0 <= p <= 1
    return LG.SimpleDiGraphFromIterator(flip_edge(rng, e, p) for e in LG.edges(graph))
end

function random_edge_orientation(
    graph::LG.SimpleGraph,
    p::Real=0.5
)
    return random_edge_orientation(Random.default_rng(), graph, p)
end

"""
    flip_edge(
        [rng::Random.AbstractRNG,]
        edge::LG.Edge,
        p::Real=0.5
    )::LG.Edge

return `edge` with probability `p` or `reverse(edge)` with probability `1-p`.
"""
function flip_edge(
    rng::Random.AbstractRNG,
    edge::LG.Edge,
    p::Real=0.5
)::LG.Edge
    return rand(rng) < p ? edge : reverse(edge)
end

flip_edge(edge::LG.Edge, p::Real=0.5) = flip_edge(Random.default_rng(), edge, p)

"""
    random_neighbor_assignment(
        [rng::Random.AbstractRNG,]
        graph::LG.SimpleGraph{T},
        roots=T[]
    )::LG.SimpleDiGraph{T} where {T}

Return a oriented subgraph of `graph` where each vertex except the `roots` is connected to a unique neighbor, i.e., each vertex has outdegree equal to one.
"""
function random_neighbor_assignment(
    rng::Random.AbstractRNG,
    graph::LG.SimpleGraph{T},
    roots=T[]
)::LG.SimpleDiGraph{T} where {T}
    g = LG.SimpleDiGraph(LG.nv(graph))
    for v in LG.vertices(g)
        if v ∉ roots
            w = rand(rng, LG.neighbors(graph, v))
            LG.add_edge!(g, v, w)
        end
    end
    return g
end

function random_neighbor_assignment(
    graph::LG.SimpleGraph{T},
    roots=T[]
) where {T}
    return random_neighbor_assignment(Random.default_rng(), graph, roots)
end

## String methods

function eachmatch_ranges(regex::Regex, string::String; overlap=false)
    matches = eachmatch(regex, string; overlap=overlap)
    p = length(regex.pattern)
    return (m.offset:(m.offset + p - 1) for m in matches)
end

function find_prefix_suffix(s::String)
    return [i for i in 1:div(length(s), 2) if s[1:i] == s[end-i+1:end]]
end
