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

window(hc::HardCorePointProcess) = hc.window
dimension(hc::HardCorePointProcess) = dimension(window(hc))

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

# Partial Rejection Sampling
# The implementation relies on Distances.pairwise, points are stored as columns of a Matrix.

"""
Sample from spatial hard core model a.k.a Poisson disk sampling with intensity λ and radius r on [0, 1]^2, using Partial Rejection Sampling [https://arxiv.org/pdf/1801.07342.pdf](H. Guo, M. Jerrum).
!!! warning "Side effects"

    In the current implementation, the return sample may contain some point outside of the box! (see generate_sample_poisson_union_disks)
"""
function generate_sample_prs(
        hc::HardCorePointProcess{T};
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
)::Vector{T} where {T}

    rng = getRNG(rng)
    window_ = win === nothing ? window(hc) : win

    n = rand(rng, Distributions.Poisson(hc.β * volume(window_)))
    points = Matrix{Float64}(undef, dimension(hc), n)
    for x in eachcol(points)
        x .= rand(window_; rng=rng)
    end

    while true
        bad = vec(any(pairwise_distances(points) .< hc.r, dims=2))
        !any(bad) && break

        resampled = generate_sample_poisson_union_disks(hc.β, hc.r, points[:, bad]; rng=rng)
        points = hcat(points[:, .!bad], resampled)
    end
    return [points[:, i] for i in 1:size(points, 2)]
end

"""
Sample from Poisson(λ) on ⋃ᵢ B(cᵢ, r) (union of disks of radius r centered at cᵢ)

Use the independence property of the Poisson point process

  - Sample from Poisson(λ) on B(c₁, r)
  - Sample from Poisson(λ) on B(c₂, r) ∖ B(c₁, r)
  - Sample from Poisson(λ) on B(cⱼ, r) ∖ ⋃_i<j B(cᵢ, r)
  - ...
"""
function generate_sample_poisson_union_disks(
        λ::Real,
        r::Real,
        centers::Matrix;
        rng=-1
)::Matrix{Float64}
    rng = getRNG(rng)
    𝒫 = Distributions.Poisson(λ * π * r^2)
    points = Matrix{Float64}(undef, 2, 0)
    for (i, c) in enumerate(eachcol(centers))
        n = rand(rng, 𝒫)
        n == 0 && continue
        proposed_points = generate_sample_uniform_in_disk(n, c, r; rng=rng)
        if i == 1
            points = hcat(points, proposed_points)
            continue
        end
        centers_ = @view centers[:, 1:i-1]
        accept = vec(all(pairwise_distances(centers_, proposed_points) .> r, dims=1))
        points = hcat(points, proposed_points[:, accept])
    end
    return points
end

"""
    generate_sample_uniform_in_disk(n, center, radius,

Sample `n` points uniformly in ``D(c, r)`` (disk with `center` and `radius`)
"""
function generate_sample_uniform_in_disk(
        n::Integer,
        center::AbstractVector,
        radius::Real;
        rng=-1
)::Matrix{Float64}
    @assert radius > 0
    @assert length(center) == 2

    n == 0 && return Matrix{Float64}(undef, 2, 0)
    rng = getRNG(rng)

    sample = repeat(center, 1, n)
    for x in eachcol(sample)
        r = radius * sqrt(rand(rng))
        theta = 2π * rand(rng)
        x .+= [r * cos(theta), r * sin(theta)]
    end
    return sample
end

# Default sample generate_sample

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
