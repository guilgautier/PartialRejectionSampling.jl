@doc raw"""
    StraussPointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}

Spatial point process with density (w.r.t. the homogenous Poisson point process with unit intensity) given by

```math
    \prod_{x \in X}
        \beta
    \prod_{\{x, y\} \subseteq X}
        \gamma^{ 1_{ \left\| x - y \right\|_2 \leq r } }
    =
    \beta^{|X|}
    \gamma^{|\{\{x, y\} \subseteq X ~;~ \left\| x - y \right\|_2 \leq r\}|}
```

with intensity ``\beta > 0``, interaction coefficient ``0\leq \gamma\leq 1`` and interaction range ``r > 0``

- ``\gamma = 0`` corresponds to the [`PRS.HardCorePointProcess`](@ref)
- ``\gamma = 1`` corresponds to the [`PRS.HomogeneousPoissonPointProcess`](@ref)

**See also**
- 6.2.2. [MoWa04](@cite)
"""
struct StraussPointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}
    "Intensity"
    β::Float64
    "Interaction coefficient"
    γ::Float64
    "Interaction range"
    r::Float64
    window::AbstractSpatialWindow{Float64}
end

@doc raw"""
    StraussPointProcess(β::Real, γ::Real, r::Real, window::AbstractSpatialWindow)

Construct a [`PRS.StraussPointProcess`](@ref) with intensity ``\beta > 0``, interaction coefficient ``0 \leq \gamma \leq 1``, and interaction range ``r > 0``, restricted to `window`.
"""
function StraussPointProcess(β::Real, γ::Real, r::Real, window::AbstractSpatialWindow)
    @assert β > 0
    @assert 0 <= γ <= 1
    @assert r > 0
    return StraussPointProcess{Vector{Float64}}(β, γ, r, window)
end

intensity(pp::StraussPointProcess) = pp.β
interaction_coefficient(pp::StraussPointProcess) = pp.γ
interaction_range(pp::StraussPointProcess) = pp.r

# Sampling

"""
    generate_sample(
        pp::StraussPointProcess{T};
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
    )::Vector{T} where {T}

Genererate an exact sample from [`PRS.StraussPointProcess`](@ref).
Default sampler is dominated CFTP of [MoWa99](@cite), see [`PRS.generate_sample_dcftp`](@ref)
"""
function generate_sample(
    pp::StraussPointProcess{T};
    win::Union{Nothing,AbstractWindow}=nothing,
    rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    return generate_sample_dcftp(pp; win=win, rng=rng)
end

## dominated CFTP, see dominated_cftp.jl

isrepulsive(pp::StraussPointProcess) = true
isattractive(pp::StraussPointProcess) = false

function papangelou_conditional_intensity(
    pp::StraussPointProcess{Vector{T}},
    x::AbstractVector{T},
    X::Union{AbstractVector{Vector{T}}, AbstractSet{Vector{T}}}
)::Real where {T}
    x in X && return 0.0
    β, γ, r = intensity(pp), interaction_coefficient(pp), interaction_range(pp)
    nb_interactions = 0
    for y in X
        nb_interactions += Distances.euclidean(x, y) <= r
    end
    return β * γ^nb_interactions
end

upper_bound_papangelou_conditional_intensity(pp::StraussPointProcess) = intensity(pp)

## Partial rejection sampling

"""
    generate_sample_prs(
        pp::StraussPointProcess{T};
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
    )::Vector{T} where {T}

Genererate an exact sample from [`PRS.StraussPointProcess`](@ref) using partial rejection sampling (PRS).
Default sampler is grid partial rejection sampling of [MoKr20](@cite), see [`PRS.generate_sample_grid_prs`](@ref)
"""
function generate_sample_prs(
    pp::StraussPointProcess{T};
    win::Union{Nothing,AbstractWindow}=nothing,
    rng=-1
)::Vector{T} where {T}
    return generate_sample_grid_prs(pp; rng=rng)
end

# grid Partial Rejection Sampling, see grid_prs.jl

@doc raw"""
    gibbs_interaction(
        pp::StraussPointProcess{T},
        cell1::SpatialCellGridPRS{T},
        cell2::SpatialCellGridPRS{T}
    )::Real where {T}

Compute the pairwise Gibbs interaction for a [`PRS.StraussPointProcess`](@ref) between ``C_1``=`cell1` and ``C_2``=`cell2`,

```math
    \prod_{(x, y) \in C_1 \times C_2}
        \gamma^{ 1_{ \left\|x - y\right\|_2 \leq r } }.
```
"""
function gibbs_interaction(
    pp::StraussPointProcess{T},
    cell1::SpatialCellGridPRS{T},
    cell2::SpatialCellGridPRS{T}
)::Real where {T}
    γ, r = interaction_coefficient(pp), interaction_range(pp)
    nb_interactions = 0
    if !isempty(cell1) && !isempty(cell2)
        for x in cell1
            for y in cell2
                nb_interactions += Distances.euclidean(x, y) <= r
            end
        end
    end
    return γ^nb_interactions
end
