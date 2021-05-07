@doc raw"""
    HardCorePointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}

Spatial point process with density (w.r.t. the homogenous Poisson point process with unit intensity) proportional to

```math
    \prod_{x \in X}
        \beta
    \prod_{\{x, y\} \subseteq X}
        1_{ \left\| x - y \right\|_2 > r },
```

where ``\beta > 0`` is called the background intensity and ``r > 0`` the interaction range.

[`PRS.HardCorePointProcess`](@ref) can be viewed as
- [`PRS.HomogeneousPoissonPointProcess`](@ref) conditioned to having no pair of points at distance less than ``r``,
- [`PRS.StraussPointProcess`](@ref) with interaction coefficient ``\gamma=0``.

**See also**

- [`PRS.HardCoreGraph`](@ref), the graph counterpart of [`PRS.HardCorePointProcess`](@ref).

# Example

A realization for ``\beta=38`` and ``r=0.1`` on ``[0, 1]^2``, where points are marked with a circle of radius ``r/2``.

![assets/hard_core_spatial.png](assets/hard_core_spatial.png)
"""
struct HardCorePointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}
    "Intensity"
    β::Float64
    "Interation range"
    r::Float64
    window::AbstractSpatialWindow{Float64}
end

function Base.show(io::IO, pp::HardCorePointProcess{T}) where {T}
    print(io, "HardCorePointProcess{$T}\n- β = $(pp.β)\n- r = $(pp.r)\n- window = $(pp.window)")
end

@doc raw"""
    HardCorePointProcess(
        β::Real,
        r::Real,
        window::Union{Nothing,AbstractSpatialWindow}=nothing
    )

Construct a [`PRS.HardCorePointProcess`](@ref) with intensity `β` and interaction range `r`, restricted to `window`.

Default window (`window=nothing`) is [`PRS.SquareWindow`](@ref)`()`.

```jldoctest; output = true
using PartialRejectionSampling

β, r = 40, 0.05
win = PRS.SquareWindow(zeros(2), 1)
PRS.HardCorePointProcess(β, r, win)

# output

HardCorePointProcess{Array{Float64,1}}
- β = 40.0
- r = 0.05
- window = SquareWindow [0.0, 1.0]^2
```
"""
function HardCorePointProcess(
    β::Real,
    r::Real,
    window::Union{Nothing,AbstractSpatialWindow}=nothing
)
    @assert β > 0
    @assert r > 0
    win = isnothing(window) ? SquareWindow() : window
    return HardCorePointProcess{Vector{Float64}}(β, r, win)
end

"""
    intensity(pp::HardCorePointProcess) = pp.β
"""
intensity(pp::HardCorePointProcess) = pp.β

"""
    interaction_coefficient(pp::HardCorePointProcess) = 0.0
"""
interaction_coefficient(pp::HardCorePointProcess) = 0.0

"""
    interaction_range(pp::HardCorePointProcess) = pp.r
"""
interaction_range(pp::HardCorePointProcess) = pp.r

## Sampling

"""
    generate_sample(
        [rng::Random.AbstractRNG,]
        pp::HardCorePointProcess{T},
        win::Union{Nothing,AbstractWindow}=nothing
    )::Vector{T} where {T}

Generate an exact sample from the [`PRS.HardCorePointProcess`](@ref).

Default window (`win=nothing`) is `window(pp)=pp.window`.

Default sampler is [`PRS.generate_sample_prs`](@ref).

**See also**

- [`PRS.generate_sample_dcftp`](@ref),
- [`PRS.generate_sample_grid_prs`](@ref).
"""
function generate_sample(
    rng::Random.AbstractRNG,
    pp::HardCorePointProcess{T},
    win::Union{Nothing,AbstractWindow}=nothing
)::Vector{T} where {T}
    return generate_sample_prs(rng, pp, win)
end

# function generate_sample(
#     pp::HardCorePointProcess,
#     win::Union{Nothing,AbstractWindow}=nothing
# )
#     return generate_sample(Random.default_rng(), pp, win)
# end

### dominated CFTP, see dominated_cftp.jl

"""
    isrepulsive(pp::HardCorePointProcess) = true
"""
isrepulsive(pp::HardCorePointProcess) = true

"""
    isattractive(pp::HardCorePointProcess) = false
"""
isattractive(pp::HardCorePointProcess) = false

@doc raw"""
    papangelou_conditional_intensity(
        pp::HardCorePointProcess{Vector{T}},
        x::AbstractVector{T},
        X::Union{AbstractVector{Vector{T}},AbstractSet{Vector{T}}}
    )::Real where {T}

Compute the [Papangelou conditional intensity](https://en.wikipedia.org/wiki/Point_process#Papangelou_intensity_function) of the point process `pp`

```math
    \beta
    \prod_{y\in X} 1_{\{\left\| x - y \right\|_2 > r\}},
```

where ``\beta=`` `pp.β` and ``r=`` `pp.r`.
"""
function papangelou_conditional_intensity(
    pp::HardCorePointProcess{Vector{T}},
    x::AbstractVector{T},
    X::Union{AbstractVector{Vector{T}},AbstractSet{Vector{T}}}
)::Real where {T}
    x in X && return 0.0
    β, r = intensity(pp), interaction_range(pp)
    return β * !any(Distances.euclidean(x, y) <= r for y in X)
end

"""
    upper_bound_papangelou_conditional_intensity(pp::HardCorePointProcess) = intensity(pp)
"""
upper_bound_papangelou_conditional_intensity(pp::HardCorePointProcess) = intensity(pp)

### Partial Rejection Sampling

@doc raw"""
    generate_sample_prs(
        [rng::Random.AbstractRNG,]
        pp::HardCorePointProcess{T},
        win::Union{Nothing,AbstractWindow}=nothing
    )::Vector{T} where {T}

Sample from [`PRS.HardCorePointProcess`](@ref) using Partial Rejection Sampling (PRS) of [GuJe18](@cite).

Default window (`win=nothing`) is `window(pp)=pp.window`.

**See also**

- [`PRS.generate_sample_dcftp`](@ref)
- [`PRS.generate_sample_grid_prs`](@ref).

# Example

A illustration of the procedure for ``\beta=38`` and ``r=0.1`` on ``[0, 1]^2`` where points are marked with a circle of radius ``r/2``.
Red disks highlight bad regions and orange disk, with radius ``r``, describe the resampling region.

![assets/hard_core_spatial_prs.gif](assets/hard_core_spatial_prs.gif)
"""
function generate_sample_prs(
    rng::Random.AbstractRNG,
    pp::HardCorePointProcess{T},
    win::Union{Nothing,AbstractWindow}=nothing
)::Vector{T} where {T}

    β, r = intensity(pp), interaction_range(pp)
    win_ = isnothing(win) ? window(pp) : win

    β_max = 0.21027 / volume(BallWindow(zeros(dimension(win_)), r/2))
    if β > β_max
        @warn "The arguments do not satisfy Lemma 6 of Guo and Jerrum (2018).\nPartial rejection sampling may not be efficient.\nGiven the interaction range `pp.r`, consider choosing β ≤ $(β_max)"
    end

    points = generate_sample(rng, HomogeneousPoissonPointProcess(β, win_))
    while true
        bad = vec(any(pairwise_distances(points) .<= r, dims=2))
        !any(bad) && break
        resampled = generate_sample_poisson_union_balls(rng, β, points[:, bad], r, win_)
        points = hcat(points[:, .!bad], resampled)
    end
    return collect(T, eachcol(points))
end

# function generate_sample_prs(
#     pp::HardCorePointProcess,
#     win::Union{Nothing,AbstractWindow}=nothing
# )
#     return generate_sample_prs(Random.default_rng(), pp, win=win)
# end

### grid Partial Rejection Sampling, see grid_prs.jl

@doc raw"""
    gibbs_interaction(
        pp::HardCorePointProcess{T},
        cell1::SpatialCellGridPRS{T},
        cell2::SpatialCellGridPRS{T}
    )::Real where {T}

Compute the pairwise Gibbs interaction for a [`PRS.HardCorePointProcess`](@ref) between ``C_1``=`cell1` and ``C_2``=`cell2`,

```math
    \prod_{(x, y) \in C_1 \times C_2}
        1_{ \left\|x - y\right\|_2 \leq r }.
```
"""
function gibbs_interaction(
    pp::HardCorePointProcess{T},
    cell1::SpatialCellGridPRS{T},
    cell2::SpatialCellGridPRS{T},
)::Real where {T}
    (isempty(cell1) || isempty(cell2)) && return 1.0
    r = interaction_range(pp)
    return !any(any(Distances.euclidean(x, y) <= r for x in cell1) for y in cell2)
end
