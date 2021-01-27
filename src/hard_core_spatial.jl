@doc raw"""
    HardCorePointProcess

``\beta > 0``, ``r > 0``,

```math
    \mathbb{P}[\mathcal{X}=X]
    \propto
    \beta^{|X|} \prod_{\{x, y\} \subseteq X} 1_{\operatorname{distance}(x, y) > r}
```

!!! seealso

    HardCorePointProcess(β, window) ≡ HardCorePointProcess(β, γ=0, window)
"""
struct HardCorePointProcess{T} <: AbstractSpatialPointProcess{T}
    β::Float64
    "Interation range"
    r::Float64
    window::AbstractSpatialWindow{Float64}
end

function HardCorePointProcess(β::Real, r::Real, window::AbstractSpatialWindow)
    @assert β > 0
    @assert r > 0
    return HardCorePointProcess{Vector{Float64}}(float(β), float(r), window)
end

function HardCorePointProcess(β, r, c, w)
    return HardCorePointProcess(β, r, spatial_window(c, w))
end

## Sampling

# Used in dominated CFTP, see dominated_cftp.jl

isrepulsive(hc::HardCorePointProcess) = true
isattractive(hc::HardCorePointProcess) = false

function papangelou_conditional_intensity(
        hc::HardCorePointProcess{Vector{T}},
        x::AbstractVector{T},
        X::Union{AbstractVector{Vector{T}},AbstractSet{Vector{T}}}
)::Real where {T}
    x in X && return 0.0
    return hc.β * !any(Distances.euclidean(x, y) <= hc.r for y in X)
end

upper_bound_papangelou_conditional_intensity(hc::HardCorePointProcess) = hc.β

# Used in grid Partial Rejection Sampling, see grid_prs.jl

function gibbs_interaction(
        hc::HardCorePointProcess{T},
        c_i::SpatialCellGridPRS{T},
        c_j::SpatialCellGridPRS{T},
)::Real where {T}
    (isempty(c_i) || isempty(c_j)) && return 1.0
    return !any(any(Distances.euclidean(x, y) <= hc.r for x in c_i) for y in c_j)
end

# Partial Rejection sampling

"""
Sample from spatial hard core model a.k.a Poisson disk sampling with intensity λ and radius r on [0, 1]^2, using Partial Rejection Sampling [https://arxiv.org/pdf/1801.07342.pdf](H. Guo, M. Jerrum).
"""
function generate_sample_prs(
        hc::HardCorePointProcess{T};
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
)::Vector{T} where {T}

    rng = getRNG(rng)
    window_ = win === nothing ? window(hc) : win

    points = generate_sample(HomogeneousPoissonPointProcess(hc.β, window_); rng=rng)

    while true
        bad = vec(any(pairwise_distances(points) .< hc.r, dims=2))
        !any(bad) && break

        resampled = generate_sample_poisson_union_balls(hc.β, points[:, bad], hc.r;
                                                        win=window_, rng=rng)
        points = hcat(points[:, .!bad], resampled)
    end
    return collect(T, eachcol(points))
end

# Default sampler generate_sample

"""
Default sampler for HardCorePointProcess
"""
function generate_sample(
        hc::HardCorePointProcess{T};
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    return generate_sample_prs(hc; win=win, rng=rng)
end
