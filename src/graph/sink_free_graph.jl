"""
    SinkFreeGraph{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}

Concrete type representing a point process on the edges of a `graph` characterizing the uniform distribution on the orientations of the edges conditioned on the absence of sinks.
"""
struct SinkFreeGraph{T<:LG.SimpleDiGraph{Int64}} <: AbstractGraphPointProcess{T}
    "Graph"
    graph::LG.SimpleGraph{Int64}
end

"""
    SinkFreeGraph(graph::LG.SimpleGraph{T}) where {T<:Int}

Construct a [`SinkFreeGraph`](@ref).
"""
SinkFreeGraph(graph::LG.SimpleGraph{Int64}) = SinkFreeGraph{LG.SimpleDiGraph{Int64}}(graph)

"""
    generate_sample(
        pp::SinkFreeGraph{T};
        rng=-1
    )::Vector{T} where {T}

Generate an exact sample from the [`PRS.SinkFreeGraph`](@ref).

Default sampler is [`PRS.generate_sample_prs`](@ref).
"""
function generate_sample(
    pp::SinkFreeGraph{T};
    rng=-1
)::Vector{T} where {T}
    return generate_sample_prs(pp; rng=rng)
end

"""
    generate_sample_prs(
        pp::SinkFreeGraph{T};
        rng=-1
    )::T where {T}

Generate an orientated version of `pp.graph` uniformly at random among all possible orientations conditioned on the absence of sinks, using Partial Rejection Sampling (PRS).

**See also**

- Section 4.1 of [GuJeLi19](@cite).
"""
function generate_sample_prs(
    pp::SinkFreeGraph{T};
    rng=-1
)::T where {T}
    rng = getRNG(rng)
    g = random_edge_orientation(pp.graph; rng=rng)
    while true
        sinks = sink_nodes(g)
        isempty(sinks) && break
        # Resample i.e. flip edges forming a sink uniformly at random (Bernoulli(0.5))
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
