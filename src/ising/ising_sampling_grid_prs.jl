## Grid PRS, see grid_prs.jl

function weighted_interaction_graph(
        ising::Ising;
        rng=-1
)::SWG.SimpleWeightedGraph
    rng = getRNG(rng)
    g = SWG.SimpleWeightedGraph(ising.g)
    for e in LG.edges(g)
        i, j = Tuple(e)
        @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng, weighttype(g))
    end
    return g
end

function initialize_cells(
        ising::Ising{T},
        size::Integer
)::Vector{GraphCellGridPRS{T}} where {T}
    return [GraphCellGridPRS(GraphNode(i), zero(T)) for i in 1:size]
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

function gibbs_interaction(
        ising::Ising{T},
        xᵢ::GraphCellGridPRS{T},
        xⱼ::GraphCellGridPRS{T}
)::Real where {T}
    return exp(ising.J * xᵢ.value * xⱼ.value - abs(ising.J))
end

function is_inner_interaction_possible(
        ising::Ising{T},
        xᵢ::GraphCellGridPRS{T},
        xⱼ::GraphCellGridPRS{T}
)::Bool where {T}
    return (sign(ising.J) * xᵢ.value * xⱼ.value) < 0
end

function is_outer_interaction_possible(
        ising::Ising{T},
        xᵢ::GraphCellGridPRS{T},
        xⱼ::GraphCellGridPRS{T}
)::Bool where {T}
    return true
end
