@doc raw"""
    StraussPointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}

Spatial point process with density (w.r.t. the homogenous Poisson point process with unit intensity) proportional to

```math
    \prod_{x \in X}
        \beta
    \prod_{\{x, y\} \subseteq X}
        \gamma^{ 1_{ \left\| x - y \right\|_2 \leq r } }
    =
    \beta^{|X|}
    \gamma^{|\{\{x, y\} \subseteq X ~;~ \left\| x - y \right\|_2 \leq r\}|},
```

with intensity ``\beta > 0``, interaction coefficient ``0\leq \gamma\leq 1`` and interaction range ``r > 0``.

- ``\gamma = 0`` corresponds to the [`PRS.HardCorePointProcess`](@ref),
- ``\gamma = 1`` corresponds to the [`PRS.HomogeneousPoissonPointProcess`](@ref).

**See also**

- Section 6.2.2 of [MoWa04](@cite).
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

function Base.show(io::IO, pp::StraussPointProcess{T}) where {T}
    print(io, "StraussPointProcess{$T}\n- β = $(pp.β)\n- γ = $(pp.γ)\n- r = $(pp.r)\n- window = $(pp.window)")
end

@doc raw"""
    StraussPointProcess(
        β::Real,
        γ::Real,
        r::Real,
        window::Union{Nothing,AbstractSpatialWindow}=nothing
    )

Construct a [`PRS.StraussPointProcess`](@ref) with intensity `β`, interaction coefficient `γ`, and interaction range `r`, restricted to `window`.

Default window (`window=nothing`) is [`PRS.SquareWindow`](@ref)`()`.

```jldoctest; output = true
using PartialRejectionSampling

β, γ, r = 2.0, 0.2, 0.7
win = PRS.SquareWindow([0.0, 0.0], 10.0)
PRS.StraussPointProcess(β, γ, r, win)

# output

StraussPointProcess{Array{Float64,1}}
- β = 2.0
- γ = 0.2
- r = 0.7
- window = SquareWindow [0.0, 10.0]^2
```

# Example

A illustration of the procedure for ``\beta=78, \gamma=0.1`` and ``r=0.07`` on ``[0, 1]^2`` where points are marked with a circle of radius ``r/2``.

![assets/strauss_spatial_prs.gif](assets/strauss_spatial_prs.gif)
"""
function StraussPointProcess(
    β::Real,
    γ::Real,
    r::Real,
    window::Union{Nothing,AbstractSpatialWindow}=nothing
)
    @assert β > 0
    @assert 0 <= γ <= 1
    @assert r > 0
    win = isnothing(window) ? SquareWindow() : window
    return StraussPointProcess{Vector{Float64}}(β, γ, r, win)
end

"""
    intensity(pp::StraussPointProcess) = pp.β
"""
intensity(pp::StraussPointProcess) = pp.β

"""
    interaction_coefficient(pp::StraussPointProcess) = pp.γ
"""
interaction_coefficient(pp::StraussPointProcess) = pp.γ

"""
    interaction_range(pp::StraussPointProcess) = pp.r
"""
interaction_range(pp::StraussPointProcess) = pp.r

# Sampling

"""
    generate_sample(
        [rng::Random.AbstractRNG,]
        pp::StraussPointProcess{T},
        win::Union{Nothing,AbstractWindow}=nothing
    )::Vector{T} where {T}

Genererate an exact sample from [`PRS.StraussPointProcess`](@ref).

Default window (`win=nothing`) is `window(pp)=pp.window`.

Default sampler is [`PRS.generate_sample_dcftp`](@ref).

**See also**
- [`PRS.generate_sample_grid_prs`](@ref).
"""
function generate_sample(
    rng::Random.AbstractRNG,
    pp::StraussPointProcess{T},
    win::Union{Nothing,AbstractWindow}=nothing
)::Vector{T} where {T}
    return generate_sample_dcftp(rng, pp, win)
end

# function generate_sample(
#     pp::StraussPointProcess,
#     win::Union{Nothing,AbstractWindow}=nothing
# )
#     return generate_sample(Random.default_rng(), pp, win)
# end

## dominated CFTP, see dominated_cftp.jl

@doc raw"""
    isrepulsive(pp::StraussPointProcess) = true

We consider only repulsive [`PRS.StraussPointProcess`](@ref) since ``0 \leq \gamma \leq 1``.
"""
isrepulsive(pp::StraussPointProcess) = true

@doc raw"""
    isattractive(pp::StraussPointProcess) = false

We consider only repulsive [`PRS.StraussPointProcess`](@ref) since ``0 \leq \gamma \leq 1``.
"""
isattractive(pp::StraussPointProcess) = false

@doc raw"""
    papangelou_conditional_intensity(
        pp::StraussPointProcess{Vector{T}},
        x::AbstractVector{T},
        X::Union{AbstractVector{Vector{T}}, AbstractSet{Vector{T}}}
    )::Real where {T}

Compute the [Papangelou conditional intensity](https://en.wikipedia.org/wiki/Point_process#Papangelou_intensity_function) of the point process `pp`

```math
    \beta
    \gamma^{|\{y \in X ~;~ \left\| x - y \right\|_2 \leq r\}|}
    ~ 1_{x \notin X},
```

where ``\beta=`` `pp.β`, ``\gamma=`` `pp.γ` and ``r=`` `pp.r`.
"""
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

"""
    upper_bound_papangelou_conditional_intensity(pp::StraussPointProcess) = intensity(pp)
"""
upper_bound_papangelou_conditional_intensity(pp::StraussPointProcess) = intensity(pp)

## Partial rejection sampling

"""
    generate_sample_prs(
        [rng::Random.AbstractRNG,]
        pp::StraussPointProcess{T},
        win::Union{Nothing,AbstractWindow}=nothing
    )::Vector{T} where {T}

Genererate an exact sample from [`PRS.StraussPointProcess`](@ref) using Partial Rejection Sampling (PRS).

Default window (`win=nothing`) is `window(pp)=pp.window`.

Default sampler is [`PRS.generate_sample_grid_prs`](@ref).
"""
function generate_sample_prs(
    rng::Random.AbstractRNG,
    pp::StraussPointProcess{T},
    win::Union{Nothing,AbstractWindow}=nothing
)::Vector{T} where {T}
    return generate_sample_grid_prs(rng, pp, win)
end

# function generate_sample_prs(
#     pp::StraussPointProcess,
#     win::Union{Nothing,AbstractWindow}=nothing
# )
#     return generate_sample_prs(Random.default_rng(), pp, win)
# end

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
