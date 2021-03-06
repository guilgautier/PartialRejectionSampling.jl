@doc raw"""
    HardCoreGraph{T<:Integer} <: AbstractGraphPointProcess{T}

Concrete type representing a point process on the vertices of a `graph` ``=(V, E)`` parametrized by ``\beta \geq 0`` which characterizes the distribution on the [independent sets](https://en.wikipedia.org/wiki/Independent_set_(graph_theory)) of `graph`, where each vertex is present with marginal probability ``\frac{\beta}{1+\beta}``.

In other words, it can also be viewed as the product distribution ``\operatorname{Bernoulli}(\frac{\beta}{1+\beta})^{\otimes |V|}`` on the vertices of `graph` conditioned on forming an independent set,

```math
    \mathbb{P}\!\left[ \mathcal{X} = X \right]
    \propto
    \prod_{x\in X}
        \frac{\beta}{1+\beta}
     1_{X \text{forms an independent set}}.
```

**See also**

- Section 7.2 of [GuJeLi19](@cite),
- Example 4.1 of [MoKr20](@cite),
- [`PRS.HardCorePointProcess`](@ref), the spatial counterpart of [`PRS.HardCoreGraph`](@ref).

# Example

A realization from a ``5\times 5`` grid graph, with ``\beta = 0.5``.

![assets/hard_core_graph.png](assets/hard_core_graph.png)
"""
struct HardCoreGraph{T<:Integer} <: AbstractGraphPointProcess{T}
    "Graph"
    graph::LG.SimpleGraph{T}
    β::Float64
end

function Base.show(io::IO, pp::HardCoreGraph{T}) where {T}
    print(io, "HardCoreGraph{$T}\n- graph = $(pp.graph)\n- β = $(pp.β)")
end

"""
    HardCoreGraph(
        graph::LG.SimpleGraph{T},
        β::Real
    ) where {T<:Integer}

Construct a [`PRS.HardCoreGraph`](@ref).

```jldoctest; output = true
using PartialRejectionSampling
using LightGraphs; const LG = LightGraphs

g, β = LG.grid([5, 5]), 1
PRS.HardCoreGraph(g, β)

# output

HardCoreGraph{Int64}
- graph = {25, 40} undirected simple Int64 graph
- β = 1.0
```
"""
function HardCoreGraph(
    graph::LG.SimpleGraph{T},
    β::Real
) where {T<:Integer}
    @assert β >= 0
    return HardCoreGraph{T}(graph, β)
end

# Sampling

"""
    generate_sample(
        [rng=Random.AbstractRNG,]
        pp::HardCoreGraph{T}
    )::Vector{T} where {T}

Generate an exact sample from the [`PRS.SinkFreeGraph`](@ref).

Default sampler is [`PRS.generate_sample_prs`](@ref).
"""
function generate_sample(
    rng::Random.AbstractRNG,
    pp::HardCoreGraph{T}
)::Vector{T} where {T}
    return generate_sample_prs(rng, pp)
end

# generate_sample(pp::HardCoreGraph) = generate_sample(Random.default_rng(), pp)

## Partial Rejeciton Sampling (PRS)

"""
    generate_sample_prs(
        [rng=Random.AbstractRNG,]
        pp::HardCoreGraph{T}
    )::Vector{T} where {T}

Sample from [`PRS.HardCoreGraph`](@ref) using Partial Rejection Sampling (PRS), see Section 7.2 of [GuJeLi19](@cite)

**See also**

- Example 4.1 of [MoKr20](@ref).

# Example

An illustration of the procedure on a ``5\\times 5`` grid graph.

- red: variables involved in constraints violation (gray neighboring nodes)
- orange: variables to be resampled (red nodes and their neighborhood)

![assets/hard_core_graph.gif](assets/hard_core_graph.gif)
"""
function generate_sample_prs(
    rng::Random.AbstractRNG,
    pp::HardCoreGraph{T}
)::Vector{T} where {T}

    β_max = 1 / (2 * sqrt(exp(1)) * LG.Δ(pp.graph) - 1)
    if pp.β > β_max
        @warn "The arguments do not satisfy Theorem 35 of Guo, Jerrum and Liu (2019).\nPartial rejection sampling may not be efficient.\nGiven `graph`, consider choosing β ≤ $(β_max)"
    end

    proba = pp.β / (one(pp.β) + pp.β)

    adj = LG.adjacency_matrix(pp.graph)
    occupied = randsubseq(rng, LG.vertices(pp.graph), proba)
    while true
        # Check if occupied vertices form an independent set
        sub_graph = LG.SimpleGraph(adj[occupied, occupied])
        LG.ne(sub_graph) == 0 && break

        independent, resample = T[], Set{T}()
        for cc in LG.connected_components(sub_graph)
            if length(cc) == 1  # Identify current independent vertices
                append!(independent, occupied[cc])
            else  # Construct the resampling set of vertices
                union!(resample, occupied[cc])
                for v in cc
                    union!(resample, LG.neighbors(pp.graph, occupied[v]))
                end
            end
        end
        randsubseq!(rng, occupied, collect(resample), proba)
        append!(occupied, independent)
    end
    return occupied
end

# generate_sample_prs(pp::HardCoreGraph) = generate_sample_prs(Random.default_rng(), pp)
