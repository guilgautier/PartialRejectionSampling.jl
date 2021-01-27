@doc raw"""
    StraussPointProcess

``\beta > 0``, ``0\leq \gamma\leq 1, r > 0``

```math
    \mathbb{P}[\mathcal{X}=X]
    \propto
    \beta^{|X|} \gamma^{|\{\{x, y\} \subseteq X ~;~ \operatorname{distance}(x, y) \leq r\}|}
```

``\gamma = 0`` corresponds to the hardcore model.
"""
struct StraussPointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}
    β::Float64
    γ::Float64
    "Interation range"
    r::Float64
    window::AbstractSpatialWindow{Float64}
end

function StraussPointProcess(β::Real, γ::Real, r::Real, window::AbstractSpatialWindow)
    @assert β > 0
    @assert 0 <= γ <= 1
    @assert r > 0
    return StraussPointProcess{Vector{Float64}}(β, γ, r, window)
end

## Sampling

# Used in dominated_cftp
isrepulsive(strauss::StraussPointProcess) = true
isattractive(strauss::StraussPointProcess) = false

function papangelou_conditional_intensity(
    strauss::StraussPointProcess{Vector{T}},
    x::AbstractVector{T},
    X::Union{AbstractVector{Vector{T}}, AbstractSet{Vector{T}}}
)::Real where {T}
    x in X && return 0.0
    nb_interactions = 0
    for y in X
        nb_interactions += Distances.euclidean(x, y) <= strauss.r
    end
    return strauss.β * strauss.γ^nb_interactions
end

upper_bound_papangelou_conditional_intensity(strauss::StraussPointProcess) = strauss.β

# Used in grid Partial Rejection Sampling, see grid_prs.jl
function gibbs_interaction(
        strauss::StraussPointProcess{T},
        cell_i::SpatialCellGridPRS{T},
        cell_j::SpatialCellGridPRS{T},
)::Real where {T}
    nb_interactions = 0
    if !isempty(cell_i) && !isempty(cell_j)
        for x in cell_i
            for y in cell_j
                nb_interactions += Distances.euclidean(x, y) < strauss.r
            end
        end
    end
    return strauss.γ^nb_interactions
end

"""
Default sampler for StraussPointProcess
"""
function generate_sample(
        strauss::StraussPointProcess{T};
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    return generate_sample_dcftp(strauss; win=win, rng=rng)
end
