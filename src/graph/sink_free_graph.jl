"""
    SinkFreeGraph{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}

Concrete type representing a point process on the edges of a `graph` characterizing the uniform distribution on the orientations of the edges conditioned on the absence of sinks.

**See also**

- Section 4.1 of [GuJeLi19](@cite).

# Example

A realization from a ``5\\times 5`` grid graph

![assets/sink_free_graph.png](assets/sink_free_graph.png)
"""
struct SinkFreeGraph{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}
    "Graph"
    graph::LG.SimpleGraph{Int64}
end

function Base.show(io::IO, pp::SinkFreeGraph{T}) where {T}
    print(io, "SinkFreeGraph{$T}\n- graph = $(pp.graph)")
end

"""
    SinkFreeGraph(graph::LG.SimpleGraph{T}) where {T<:Int}

Construct a [`PRS.SinkFreeGraph`](@ref).

```jldoctest; output = true
using PartialRejectionSampling
using LightGraphs; const LG = LightGraphs

sfg = PRS.SinkFreeGraph(LG.grid([5, 5]))

# output

SinkFreeGraph{LightGraphs.SimpleGraphs.SimpleDiGraph{Int64}}
- graph = {25, 40} undirected simple Int64 graph
```
"""
SinkFreeGraph(graph::LG.SimpleGraph{Int64}) = SinkFreeGraph{LG.SimpleDiGraph{Int64}}(graph)

"""
    generate_sample([rng::Random.AbstractRNG,] pp::SinkFreeGraph)

Generate an exact sample from the [`PRS.SinkFreeGraph`](@ref).

Default sampler is [`PRS.generate_sample_prs`](@ref).
"""
generate_sample(rng::Random.AbstractRNG, pp::SinkFreeGraph) = generate_sample_prs(rng, pp)
# generate_sample(pp::SinkFreeGraph) = generate_sample_prs(Random.default_rng(), pp)

"""
    generate_sample_prs(
        [rng::Random.AbstractRNG,]
        pp::SinkFreeGraph{T}
    )::T where {T}

Generate an orientated version of `pp.graph` uniformly at random among all possible orientations conditioned on the absence of sinks, using Partial Rejection Sampling (PRS).

**See also**

- Section 4.1 of [GuJeLi19](@cite).

# Example

A illustration of the procedure on a ``5 \\times 5`` grid grah.

- red: variables involved in constraints violation (orientation of the edges pointing to a sink node)
- orange: variables resampled (orientation of the previous edges)

![assets/sink_free_graph_prs.gif](assets/sink_free_graph_prs.gif)
"""
function generate_sample_prs(
    rng::Random.AbstractRNG,
    pp::SinkFreeGraph{T}
)::T where {T}
    g = random_edge_orientation(rng, pp.graph)
    while true
        sinks = sink_nodes(g)
        isempty(sinks) && break
        # Resample i.e. flip edges forming a sink uniformly at random (Bernoulli(0.5))
        for v in sinks
            in_neighbors_v = copy(LG.inneighbors(g, v))
            for w in in_neighbors_v
                if rand(rng) < 0.5
                    LG.rem_edge!(g, w, v)
                    LG.add_edge!(g, v, w)
                end
            end
        end
    end
    return g
end

# generate_sample_prs(pp::SinkFreeGraph) = generate_sample_prs(Random.default_rng(), pp)
