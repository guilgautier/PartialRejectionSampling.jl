# Grid PRS, see grid_prs.jl

"""
    weighted_interaction_graph(
        ising::Ising;
        rng=-1
    )::SWG.SimpleWeightedGraph

Return a weighted version of `ising.graph` where each edge is attached an independent uniform random variable.

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
"""
function weighted_interaction_graph(
    ising::Ising;
    rng=-1
)::SWG.SimpleWeightedGraph
    return uniform_weighted_graph(ising.graph; rng=rng)
end

@doc raw"""
    initialize_cells(
        ising::Ising{T}
    )::Vector{GraphCellGridPRS{T}} where {T}

Each node of `ising.graph` is considered as a `cell` of type [`PRS.GraphCellGridPRS`](@ref) such that `cell.window` is a [`PRS.GraphNode`](@ref) and `cell.value` initialized to `zero(T)`.

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
"""
function initialize_cells(
    ising::Ising{T}
)::Vector{GraphCellGridPRS{T}} where {T}
    return [GraphCellGridPRS(GraphNode(i), zero(T)) for i in 1:LG.nv(ising.graph)]
end

"""
    generate_sample(
        pp::Ising;
        win::GraphNode,
        rng=-1
    )

Generate an exact sample from the marginal distribution of `pp` at state indexed by `win.idx`.
"""
function generate_sample(
    pp::Ising;
    win::GraphNode,
    rng=-1
)
    return generate_sample(pp, win.idx; rng=rng)
end

@doc raw"""
    gibbs_interaction(
        ising::Ising{T},
        xᵢ::GraphCellGridPRS{T},
        xⱼ::GraphCellGridPRS{T}
    )::Real where {T}

Compute the Gibbs interaction of [`PRS.Ising`](@ref)

```math
    \exp(J x_i x_j - |J|)).
```

**Note** the Gibbs interaction is normalized in be in [0, 1] to fit the framework of [MoKr20](@cite) and [FeViYi19](@cite).
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

Assume `xᵢ` and `xⱼ` are neighboring sites in the weighted interaction graph constructed by [`PRS.weighted_interaction_graph`](@ref) from `ising.graph` and already identified in the set of variables to be resampled [`PRS.generate_sample_grid_prs`](@ref).
Given the states of `xᵢ` and `xⱼ`, check whether a new assigment of ``U_{ij}`` (the weight of edge ``\{i,j\}``) can induce the bad event

```math
    U_{ij} > \exp(J x_i x_j - |J|))
```

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
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

Assume `xᵢ` and `xⱼ` are neighboring sites in the weighted interaction graph constructed by [`PRS.weighted_interaction_graph`](@ref) from `ising.graph`.
Given the state of `xᵢ`, one can always find a new assigment of `xⱼ` and/or ``U_{ij}`` (the weight of edge ``\{i,j\}``) can induce the bad event

```math
    U_{ij} > \exp(J x_i x_j - |J|))
```

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
"""
function is_outer_interaction_possible(
    ising::Ising{T},
    xᵢ::GraphCellGridPRS{T},
    xⱼ::GraphCellGridPRS{T}
)::Bool where {T}
    return true
end
