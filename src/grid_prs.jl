"""
Implementation of Grid Partial Rejection Sampling of [Moka, Sarat B. and Kroese, Dirk P. (2020)](https://espace.library.uq.edu.au/view/UQ:d924abb)
"""

abstract type AbstractCellGridPRS end

window(cell::AbstractCellGridPRS) = cell.window
dimension(cell::AbstractCellGridPRS) = dimension(window(cell))

mutable struct GraphCellGridPRS{T} <: AbstractCellGridPRS
    window::GraphNode{T}
    value::T
end

mutable struct SpatialCellGridPRS{T<:Vector{Float64}} <: AbstractCellGridPRS
    window::Union{RectangleWindow,SquareWindow}
    value::Vector{T}
end

Base.isempty(cell::SpatialCellGridPRS) = isempty(cell.value)
Base.iterate(cell::AbstractCellGridPRS, state=1) = iterate(cell.value, state)

function generate_sample_grid_prs(
        pp::AbstractPointProcess{T};
        rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    g = weighted_interaction_graph(pp; rng=rng)

    cells = initialize_cells(pp, LG.nv(g))
    resample_indices = Set(1:length(cells))

    while !isempty(resample_indices)
        generate_sample!(cells, resample_indices, pp; rng=rng)
        resample_indices = find_cells_to_resample_indices!(g, cells, pp; rng=rng)
    end
    return vcat(getfield.(cells, :value)...)
end

@doc raw"""
The dependency graph between the cells is a weighted king graph
where each edge ``\{i,j\}`` is associated to an event and carries a uniform mark ``U_{ij}``.
Note https://en.wikipedia.org/wiki/King%27s_graph
"""
function weighted_interaction_graph(
        pp::AbstractSpatialPointProcess;
        rng=-1
)::SWG.SimpleWeightedGraph
    rng = getRNG(rng)
    g = SWG.SimpleWeightedGraph(king_graph(ceil(Int, inv(pp.r))))
    for e in LG.edges(g)
        i, j = Tuple(e)
        @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng, weighttype(g))
    end
    return g
end

function generate_sample!(
        cell::AbstractCellGridPRS,
        pp::AbstractPointProcess;
        rng=-1
)
    rng = getRNG(rng)
    cell.value = generate_sample(pp; win=cell.window, rng=rng)
end

function generate_sample!(
        cells::Vector{T},
        indices,
        pp::AbstractPointProcess;
        rng=-1
) where {T<:AbstractCellGridPRS}
    rng = getRNG(rng)
    for i in indices
        generate_sample!(cells[i], pp; rng=rng)
    end
end

@doc raw"""
Identify bad events and return the corresponding cells' index.
An event ``\{i,j\}`` is said to be \"bad\"

```math
    \left\{U_{ij} > \exp \left[ -\sum_{x \in C_i} \sum_{y \in C_j} V(x,y) \right] \right\}
```

where ``U_{ij}`` is stored as the weight of edge ``\{i,j\}`` in the dependency graph ``g``.
Note: when a bad event occurs, the corresponding ``U_{ij}`` is resampled hence the "!"
"""
function find_bad_cells_indices!(
        g::SWG.SimpleWeightedGraph{T,U},
        cells::Vector{V},
        pp::AbstractPointProcess;
        rng=-1
)::Set{T} where {T,U,V<:AbstractCellGridPRS}
    rng = getRNG(rng)
    bad = Set{T}()
    for e ∈ LG.edges(g)
        i, j = Tuple(e)
        if e.weight > gibbs_interaction(pp, cells[i], cells[j])
            union!(bad, i, j)
            # resample U_ij associated to the bad event
            @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng, U)
        end
    end
    return bad
end

"""
Identify which events need to be resampled and return the corresponding cells' index.
This is the core of the Partial Rejection Sampling algorithm.
Note: when an event needs to be resampled, the corresponding mark ``U_{ij}`` is resampled hence the "!"
"""
function find_cells_to_resample_indices!(
        g::SWG.SimpleWeightedGraph{T,U},
        cells::Vector{V},
        pp::AbstractPointProcess;
        rng=-1
)::Set{T} where {T,U,V<:AbstractCellGridPRS}
    rng = getRNG(rng)
    R = find_bad_cells_indices!(g, cells, pp; rng=rng)
    ∂R, ∂R_tmp = copy(R), empty(R)
    while !isempty(∂R)
        for i ∈ ∂R
            isempty(cells[i]) && continue
            for j ∈ LG.neighbors(g, i)
                if j ∈ R
                    if is_inner_interaction_possible(pp, cells[i], cells[j])
                        @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng, U)
                    end
                elseif is_outer_interaction_possible(pp, cells[i], cells[j])
                    union!(R, j)
                    union!(∂R_tmp, j)
                    @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng, U)
                end
            end
        end
        ∂R, ∂R_tmp = ∂R_tmp, empty(∂R)
    end
    return R
end

## Spatial point processes

function initialize_cells(
        spp::AbstractSpatialPointProcess{T},
        size::Integer
)::Vector{SpatialCellGridPRS{T}} where {T}
    cells = Vector{SpatialCellGridPRS{T}}(undef, size)
    k = ceil(Int, inv(spp.r))
    for i in eachindex(cells)
        c_y, c_x = divrem(i-1, k)
        c = spp.window.c + spp.r .* [c_x, c_y]
        win_i = rectangle_square_window(c, min.(spp.r, spp.window.w .- c))
        cells[i] = SpatialCellGridPRS(win_i, T[])
    end
    return cells
end

function is_inner_interaction_possible(
        spp::AbstractSpatialPointProcess,
        cell_i::SpatialCellGridPRS,
        cell_j::SpatialCellGridPRS
)::Bool
    return false
end

@doc raw"""
Assuming ``C_i`` and ``C_j`` are neighboring cells,
Given the configuration of ``C_i``, check whether an assignment of ``C_j`` can induce a bad event ``\{i, j\}``.
"""
function is_outer_interaction_possible(
        spp::AbstractSpatialPointProcess,
        cell_i::SpatialCellGridPRS,
        cell_j::SpatialCellGridPRS
)::Bool
    isempty(cell_i) && return false

    win_i, win_j = window(cell_i), window(cell_j)
    # If i-j is a vertical or horizontal neighborhood
    any(win_i.c .== win_j.c) && return true

    # If i-j is a diagonal neighborhood
    c_ij = any(win_i.c .< win_j.c) ? win_j.c : win_i.c
    return any(Distances.euclidean(c_ij, x) < spp.r for x in cell_i)
end
