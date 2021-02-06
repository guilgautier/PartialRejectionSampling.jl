@doc raw"""
    HardCoreGraph{T<:Integer} <: AbstractGraphPointProcess{T}

Point process defined on the vertices of a `graph` ``=(V, E)``, and parametrized by β ``\geq 0``, with joint density

```math
    \mathbb{P}\!\left[ \mathcal{X} = X \right]
    \propto
    \prod_{x\in X}
        \frac{\beta}{1+\beta}
     1_{X \text{forms an independent set}}
```

It can be viewed as a distribution on the [independent sets](https://en.wikipedia.org/wiki/Independent_set_(graph_theory)) of `graph`, where each vertex is present with marginal probability ``\frac{\beta}{1+\beta}``.
In other words, it can also be viewed as the product distribution ``\operatorname{Bernoulli}(\frac{\beta}{1+\beta})^{\otimes |V|}`` on the vertices of `graph` conditioned on forming an independent set.

**See also**

- Section 7.2 of [GuJeLi19](@cite)
- Example 4.1 of [MoKr20](@ref)
- [`PRS.HardCoreSpatial`](@ref)
"""
struct HardCoreGraph{T<:Integer} <: AbstractGraphPointProcess{T}
    "Graph"
    graph::LG.SimpleGraph{T}
    β::Float64
end

"""
    HardCoreGraph(
        graph::LG.SimpleGraph{T},
        β::Real
    ) where {T<:Integer}

Construct a [`HardCoreGraph`](@ref) model on `graph` with parameter `β≥0`
"""
function HardCoreGraph(
    graph::LG.SimpleGraph{T},
    β::Real
) where {T<:Integer}
    @assert β >= 0
    return HardCoreGraph{T}(graph, float(β))
end

# Sampling

## Partial Rejeciton Sampling (PRS)

"""
    generate_sample_prs(
            pp::HardCoreGraph{T};
            win::Union{Nothing,AbstractWindow}=nothing,
            rng=-1
    )::Vector{T} where {T}

Sample from [`PRS.HardCoreGraph`](@ref) using Partial Rejection Sampling (PRS), see Section 7.2 of [GuJeLi19](@cite)

**See also**

- Example 4.1 of [MoKr20](@ref)
"""
function generate_sample_prs(
    hcg::HardCoreGraph{T};
    rng=-1
)::Vector{T} where {T}

    proba = hcg.β / (one(hcg.β) + hcg.β)

    rng = getRNG(rng)
    adj = LG.adjacency_matrix(hcg.graph)
    occupied = randsubseq(rng, LG.vertices(hcg.graph), proba)
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
                    union!(resample, LG.neighbors(hcg.graph, occupied[v]))
                end
            end
        end
        randsubseq!(rng, occupied, collect(resample), proba)
        append!(occupied, independent)
    end
    return occupied
end
