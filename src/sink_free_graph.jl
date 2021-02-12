"""
    SinkFreeGraph{T<:LG.SimpleGraphs.SimpleEdge{Int64}} <: AbstractGraphPointProcess{T}
"""
struct SinkFreeGraph{T<:LG.SimpleGraphs.SimpleEdge{Int64}} <: AbstractGraphPointProcess{T}
    "Graph"
    graph::LG.SimpleGraph{Int64}
end

"""
    SinkFreeGraph(graph::LG.SimpleGraph{T}) where {T<:Int}

Construct a [`SinkFreeGraph`](@ref)
"""
function SinkFreeGraph(
    graph::LG.SimpleGraph{T}
) where {T<:Int}
    return SinkFreeGraph{LG.SimpleGraphs.SimpleEdge{T}}(graph)
end

"""
    generate_sample_prs(
        sfg::SinkFreeGraph{T};
        rng=-1
    )::Vector{T} where {T}

Generate an orientated version of `sfg.graph` uniformly at random among all possible orientation where no sink occurs, using Partial Rejection Sampling (PRS), see Section 4.1 of [GuJeLi19](@cite)
"""
function generate_sample_prs(
    sfg::SinkFreeGraph{T};
    rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    g = random_edge_orientation(sfg.graph; rng=rng)
    while true
        sinks = sink_nodes(g)
        isempty(sinks) && break
        # Resample i.e. flip edges forming a sink uniformly at random
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
