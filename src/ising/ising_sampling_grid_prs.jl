# Grid PRS, see grid_prs.jl

"""
    weighted_interaction_graph(
        ising::Ising;
        rng=-1
    )::SWG.SimpleWeightedGraph

Return a weighted version of `ising.graph` where the weight of each edge is a uniform random variable.
"""
function weighted_interaction_graph(
    ising::Ising;
    rng=-1
)::SWG.SimpleWeightedGraph
    rng = getRNG(rng)
    g = SWG.SimpleWeightedGraph(ising.graph)
    for e in LG.edges(g)
        i, j = Tuple(e)
        @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng, weighttype(g))
    end
    return g
end

@doc raw"""
    initialize_cells(
        ising::Ising{T}
    )::Vector{GraphCellGridPRS{T}} where {T}

Each node of `ising.graph` is considered as `cell` of type [`GraphCellGridPRS`](@ref) such that `cell.window` is a [`GraphNode`](@ref) and `cell.value` initialized to `zero(T)`

This function is used as a subroutine of [`generate_sample_grid_prs`](@ref)
"""
function initialize_cells(
    ising::Ising{T}
)::Vector{GraphCellGridPRS{T}} where {T}
    return [GraphCellGridPRS(GraphNode(i), zero(T)) for i in 1:LG.nv(ising.graph)]
end

function generate_sample!(
    cell::GraphCellGridPRS,
    ising::Ising;
    rng=-1
)
    rng = getRNG(rng)
    hᵢ = ising.h isa Number ? ising.h : ising.h[cell.win.idx]
    cell.value = rand(rng) < sigmoid(hᵢ) ? 1 : -1
end

@doc raw"""
    gibbs_interaction(
        ising::Ising{T},
        xᵢ::GraphCellGridPRS{T},
        xⱼ::GraphCellGridPRS{T}
    )::Real where {T}

Compute the Gibbs interaction of [`Ising`](@ref)
```math
    \exp(J x_i x_j - |J|))

```

**Note**

The Gibbs interaction is normalized in be in [0, 1] to fit the framework of [MoKr20](@cite) and [FeViYi](@cite)
"""
function gibbs_interaction(
    ising::Ising{T},
    xᵢ::GraphCellGridPRS{T},
    xⱼ::GraphCellGridPRS{T}
)::Real where {T}
    return exp(ising.J * xᵢ.value * xⱼ.value - abs(ising.J))
end

@doc raw"""
    is_inner_interaction_possible(
        ising::Ising{T},
        xᵢ::GraphCellGridPRS{T},
        xⱼ::GraphCellGridPRS{T}
    ) = true

Assume `xᵢ` and `xⱼ` are neighboring sites in the weighted interaction graph constructed by [`weighted_interaction_graph`](@ref) from `ising.graph` and already identified in the set of variables to be resampled [`generate_sample_grid_prs`](@ref).
Given the states of `xᵢ` and `xⱼ`, check whether a new assigment of ``U_{ij}`` (the weight of edge ``\{i,j\}``) can induce the bad event

```math
    U_{ij} > \exp(J x_i x_j - |J|))
```

This function is used as a subroutine of [`generate_sample_grid_prs`](@ref)
"""
function is_inner_interaction_possible(
    ising::Ising{T},
    xᵢ::GraphCellGridPRS{T},
    xⱼ::GraphCellGridPRS{T}
)::Bool where {T}
    return (sign(ising.J) * xᵢ.value * xⱼ.value) < 0
end

@doc raw"""
    is_outer_interaction_possible(
        ising::Ising{T},
        xᵢ::GraphCellGridPRS{T},
        xⱼ::GraphCellGridPRS{T}
    ) = true

Assume `xᵢ` and `xⱼ` are neighboring sites in the weighted interaction graph constructed by [`weighted_interaction_graph`](@ref) from `ising.graph`.
Given the state of `xᵢ`, one can always find a new assigment of `xⱼ` and/or ``U_{ij}`` (the weight of edge ``\{i,j\}``) can induce the bad event

```math
    U_{ij} > \exp(J x_i x_j - |J|))
```

This function is used as a subroutine of [`generate_sample_grid_prs`](@ref)
"""
function is_outer_interaction_possible(
    ising::Ising{T},
    xᵢ::GraphCellGridPRS{T},
    xⱼ::GraphCellGridPRS{T}
)::Bool where {T}
    return true
end
