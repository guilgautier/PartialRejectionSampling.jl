"""
Implementation of Grid Partial Rejection Sampling of [MoKr20](@cite).
"""

"""
    AbstractCellGridPRS

Abstract type describing a cell used Grid Partial Rejection Sampling of [MoKr20](@cite).
"""
abstract type AbstractCellGridPRS end

"""
    GraphCellGridPRS{T} <: AbstractCellGridPRS

Mutable struct describing a graph node with fields

- `window` of type [`PRS.GraphNode`](@ref)`{T}`
- `value` of type `T`
"""
mutable struct GraphCellGridPRS{T} <: AbstractCellGridPRS
    window::GraphNode{T}
    value::T
end

"""
    SpatialCellGridPRS{T<:Vector{Float64}} <: AbstractCellGridPRS

Mutable struct describing a spatial cell with fields

- `window` of type `Union{`[`PRS.RectangleWindow`](@ref)`,`[`PRS.RectangleWindow`](@ref)`}`
- `value` of type `Vector{T}`
"""
mutable struct SpatialCellGridPRS{T<:Vector{Float64}} <: AbstractCellGridPRS
    window::Union{RectangleWindow,SquareWindow}
    value::Vector{T}
end

"""
    window(cell::AbstractCellGridPRS) = cell.window
"""
window(cell::AbstractCellGridPRS) = cell.window

"""
    dimension(cell::AbstractCellGridPRS) = dimension(window(cell))
"""
dimension(cell::AbstractCellGridPRS) = dimension(window(cell))

"""
    Base.iterate(cell::AbstractCellGridPRS, state=1) = iterate(cell.value, state)
"""
Base.iterate(cell::AbstractCellGridPRS, state=1) = iterate(cell.value, state)

"""
    Base.isempty(cell::SpatialCellGridPRS) = isempty(cell.value)
"""
Base.isempty(cell::SpatialCellGridPRS) = isempty(cell.value)

"""
    generate_sample!(
        rng::Random.AbstractRNG,
        cell::AbstractCellGridPRS,
        pp::AbstractPointProcess
    )

Generate an exact sample of `pp` in `cell.window` and save it in `cell.value`
"""
function generate_sample!(
    rng::Random.AbstractRNG,
    cell::AbstractCellGridPRS,
    pp::AbstractPointProcess
)
    cell.value = generate_sample(rng, pp, cell.window)
end

"""
    generate_sample!(
        rng::Random.AbstractRNG,
        cells::Vector{T},
        indices,
        pp::AbstractPointProcess
    ) where {T<:AbstractCellGridPRS}

Apply [`PRS.generate_sample!`](@ref) to each cell of `cells` indexed by `indices`.
"""
function generate_sample!(
    rng::Random.AbstractRNG,
    cells::Vector{T},
    indices,
    pp::AbstractPointProcess
) where {T<:AbstractCellGridPRS}
    for i in indices
        generate_sample!(rng, cells[i], pp)
    end
end

"""
    generate_sample_grid_prs(
        [rng::Random.AbstractRNG,]
        pp::AbstractPointProcess{T}
    )::Vector{T} where {T}

Generate an exact sample from `pp` using grid Partial Rejection Sampling (grid PRS) of [MoKr20](@cite).
"""
function generate_sample_grid_prs(
    rng::Random.AbstractRNG,
    pp::AbstractPointProcess{T}
)::Vector{T} where {T}
    g = weighted_interaction_graph(rng, pp)
    cells = initialize_cells(pp)
    resample_indices = Set(1:length(cells))

    while !isempty(resample_indices)
        generate_sample!(rng, cells, resample_indices, pp)
        resample_indices = find_cells_to_resample_indices!(rng, g, cells, pp)
    end
    return vcat(getfield.(cells, :value)...)
end

function generate_sample_grid_prs(pp::AbstractPointProcess, args...)
    return generate_sample_grid_prs(Random.default_rng(), pp, args...)
end

@doc raw"""
    gibbs_interaction(
        pp::AbstractPointProcess{T},
        cell1::AbstractCellGridPRS{T},
        cell2::AbstractCellGridPRS{T}
    )::Real where {T}

Compute the pairwise Gibbs interaction `pp` between ``C_1``=`cell1` and ``C_2``=`cell2`.

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
"""
function gibbs_interaction end

@doc raw"""
    find_bad_cells_indices!(
        [rng::Random.AbstractRNG,]
        g::SWG.SimpleWeightedGraph{T,U},
        cells::Vector{V},
        pp::AbstractPointProcess
    )::Set{T} where {T,U,V<:AbstractCellGridPRS}

Identify bad events and return the corresponding cells' index.
An event ``\{i,j\}`` is said to be \"bad\"

```math
    \left\{U_{ij} > \exp \left[ -\sum_{x \in C_i} \sum_{y \in C_j} V(x,y) \right] \right\}
```

where ``U_{ij}`` is the weight of the edge ``\{i,j\}`` in the interaction graph ``g`` created by [`PRS.weighted_interaction_graph`](@ref) and ``V`` the Gibbs potential discribing the pairwise Gibbs interaction of `pp`

**Note** when a bad event occurs, the corresponding ``U_{ij}`` is resampled hence the "!"

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
"""
function find_bad_cells_indices!(
    rng::Random.AbstractRNG,
    g::SWG.SimpleWeightedGraph{T,U},
    cells::Vector{V},
    pp::AbstractPointProcess
)::Set{T} where {T,U,V<:AbstractCellGridPRS}
    bad = Set{T}()
    for e in LG.edges(g)
        i, j = Tuple(e)
        if e.weight > gibbs_interaction(pp, cells[i], cells[j])
            union!(bad, i, j)
            # resample U_ij associated to the bad event
            @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng, U)
        end
    end
    return bad
end

@doc raw"""
    find_cells_to_resample_indices!(
        rng::Random.AbstractRNG,
        g::SWG.SimpleWeightedGraph{T,U},
        cells::Vector{V},
        pp::AbstractPointProcess
    )::Set{T} where {T,U,V<:AbstractCellGridPRS}

Identify the set of events to be resampled as constructed by Algorithm 5 in [GuJeLi19](@cite) as part of the Partial Rejection Sampling (PRS) method.
Return the indices of the variables (here cells) involved in the corresponding events.

This function is used as a subroutine of the grid PRS methodology of [MoKr20](@cite), see [`PRS.generate_sample_grid_prs`](@ref).

**Note** if the event associated to the edge ``\{i,j\}`` of `g` is selected to be resampled, the uniform random variable encoded as the weight of the correspond edge is resampled (hence the "!")

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
"""
function find_cells_to_resample_indices!(
    rng::Random.AbstractRNG,
    g::SWG.SimpleWeightedGraph{T,U},
    cells::Vector{V},
    pp::AbstractPointProcess
)::Set{T} where {T,U,V<:AbstractCellGridPRS}
    R = find_bad_cells_indices!(rng, g, cells, pp)
    ∂R, ∂R_tmp = copy(R), empty(R)
    while !isempty(∂R)
        for i in ∂R
            isempty(cells[i]) && continue
            for j in LG.neighbors(g, i)
                if j in R
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

"""
    weighted_interaction_graph(
        rng::Random.AbstractRNG,
        pp::AbstractSpatialPointProcess;
    )::SWG.SimpleWeightedGraph

Construct the weighted interaction graph ([King graph](https://en.wikipedia.org/wiki/King%27s_graph)) used in [`PRS.generate_sample_grid_prs`](@ref), to generate exact samples from [`PRS.AbstractSpatialPointProcess`](@ref)

The `pp.window` is divided into cells of length the interaction range `pp.r`.
Each cell represents a vertex of the interaction (king) graph and each edge carries a uniform random varialble.

**See also**

- Figure 4 [MoKr20](@cite)
"""
function weighted_interaction_graph(
    rng::Random.AbstractRNG,
    spp::AbstractSpatialPointProcess
)::SWG.SimpleWeightedGraph
    window_ = window(spp)
    if allequal(window_.w)
        nb_cells_x = ceil(Int, window_.w[1] / spp.r)
        return uniform_weighted_graph(rng, king_graph(nb_cells_x))
    else
        nb_cells_x, nb_cells_y = ceil.(Int, window_.w ./ spp.r)
        return uniform_weighted_graph(rng, king_graph(nb_cells_x, nb_cells_y))
    end
end

@doc raw"""
    initialize_cells(
        spp::AbstractSpatialPointProcess{T},
    )::Vector{SpatialCellGridPRS{T}} where {T}

The `spp.window` ([`PRS.RectangleWindow`](@ref) or [`PRS.SquareWindow`](@ref)) is divided into [`PRS.SpatialCellGridPRS`](@ref) of length the interaction range `pp.r` following the construction of [MoKr20](@cite) in their grid Partial Rejection Sampling (grid PRS) methodology.

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
"""
function initialize_cells(
    spp::AbstractSpatialPointProcess{T},
)::Vector{SpatialCellGridPRS{T}} where {T}
    dimension(spp) != 2 && throw(DomainError(spp, "must be 2D point process for now"))
    window_ = window(spp)
    nb_cells_x = nb_cells_y = 1
    if allequal(window_.w)  # SquareWindow
        nb_cells_x = nb_cells_y = ceil(Int, window_.w[1] / spp.r)
    else
        nb_cells_x, nb_cells_y = ceil.(Int, window_.w ./ spp.r)
    end

    cells = Vector{SpatialCellGridPRS{T}}(undef, nb_cells_x * nb_cells_y)
    for i in eachindex(cells)
        c_ij = divrem(i-1, nb_cells_x)  # coordinate (i, j) of the cell in the grid
        c_xy = @. window_.c + spp.r * c_ij
        win_i = rectangle_square_window(c_xy, @. min(spp.r, window_.w - c_xy))
        cells[i] = SpatialCellGridPRS(win_i, T[])
    end
    return cells
end

@doc raw"""
    is_inner_interaction_possible(
        spp::AbstractSpatialPointProcess,
        cell_i::SpatialCellGridPRS,
        cell_j::SpatialCellGridPRS
    ) = false

Assume `cell_i` and `cell_j` are neighboring cells in the weighted interaction graph constructed by [`PRS.weighted_interaction_graph`](@ref) from `spp` and already identified in the set of variables to be resampled in [`PRS.generate_sample_grid_prs`](@ref)
Since the configuration of points in the corresponding cells and the uniform random variable associated to the event `\{i,j\}` are considered fixed, there is no degree of freedom to make the interaction between `cell_i` and `cell_j` possible.

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
"""
function is_inner_interaction_possible(
        spp::AbstractSpatialPointProcess,
        cell_i::SpatialCellGridPRS,
        cell_j::SpatialCellGridPRS
)::Bool
    return false
end

@doc raw"""
    is_outer_interaction_possible(
        spp::AbstractSpatialPointProcess,
        cell_i::SpatialCellGridPRS,
        cell_j::SpatialCellGridPRS
    )::Bool

Assume `cell_i` and `cell_j` are neighboring cells in the weighted interaction graph constructed by [`PRS.weighted_interaction_graph`](@ref) from `spp`.
Given the configuration of points in `cell_i`, check whether a realization of `spp` in `cell_j` can induce a bad event.

This is a subroutine of [`PRS.generate_sample_grid_prs`](@ref).
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
