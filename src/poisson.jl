@doc raw"""
    HomogeneousPoissonPointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}

Homegeneous [Poisson point process](https://en.wikipedia.org/wiki/Poisson_point_process) with intensity ``\beta``, denoted ``operatorname{Poisson}(\beta)``.

``operatorname{Poisson}(\beta)`` has density (w.r.t. the homogenous Poisson point process with unit intensity ``operatorname{Poisson}(1)``) proportional to

```math
    \prod_{x \in X}
        \beta
````
"""
struct HomogeneousPoissonPointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}
    β::Float64
    window::AbstractSpatialWindow{Float64}
end

@doc raw"""
    HomogeneousPoissonPointProcess(β::Real, window::AbstractSpatialWindow)

Construct a [`PRS.HomogeneousPoissonPointProcess`](@ref) ``\operatorname{Poisson}(\beta)`` with intensity ``\beta`` restricted to `window`
"""
function HomogeneousPoissonPointProcess(β::Real, window::AbstractSpatialWindow)
    @assert β > 0
    return HomogeneousPoissonPointProcess{Vector{Float64}}(β, window)
end

intensity(hp::HomogeneousPoissonPointProcess) = hp.β

## Sampling

"""
    generate_sample(
        hp::HomogeneousPoissonPointProcess;
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
    )::Matrix{Float64}

Generate a sample from a homogenous [`PRS.HomogeneousPoissonPointProcess`](@ref) `hp` on window `win`.
Default window (`win=nothing`) is `hp.window` otherwise `win=nothing`.
"""
function generate_sample(
        hp::HomogeneousPoissonPointProcess;
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
)::Matrix{Float64}
    rng = getRNG(rng)
    window_ = win === nothing ? window(hp) : win
    n = rand(rng, Distributions.Poisson(hp.β * volume(window_)))
    return rand(window_, n; rng=rng)
end

@doc raw"""
    generate_sample_poisson_union_balls(
        β::Real,
        centers::Matrix,
        radius::Real;
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
    )::Matrix{Float64}

Generate a sample from a homogenous [`PRS.HomogeneousPoissonPointProcess`](@ref) Poisson(β) on ``\bigcup_{i} B(c_i, r)`` (union of balls centered at ``c_i`` with the same radius ``r``)

Use the independence property of the Poisson point process on disjoint subsets in order to

  - Sample from Poisson(β) on ``B(c_1, r)``
  - Sample from Poisson(β) on ``B(c_2, r) \setminus B(c_1, r)``
  - Sample from Poisson(β) on ``B(c_j, r) \setminus \bigcup_{i<j} B(c_i, r)``
  - ...
"""
function generate_sample_poisson_union_balls(
        β::Real,
        centers::Matrix,
        radius::Real;
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
)::Matrix{Float64}
    rng = getRNG(rng)
    d = size(centers, 1)
    !(win === nothing) && @assert dimension(win) == d

    𝒫 = Distributions.Poisson(β * volume(BallWindow(zeros(d), radius)))
    points = Matrix{Float64}(undef, d, 0)
    for (i, c) in enumerate(eachcol(centers))
        n = rand(rng, 𝒫)
        n == 0 && continue
        proposed = rand(BallWindow(c, radius), n; rng=rng)
        if i == 1
            points = hcat(points, proposed)
            continue
        end
        centers_ = @view centers[:, 1:i-1]
        accept = vec(all(pairwise_distances(centers_, proposed) .> radius, dims=1))
        points = hcat(points, proposed[:, accept])
    end

    if win === nothing || isempty(points)
        return points
    else
        return points[:, vec(mapslices(x -> x in win, points; dims=1))]
    end
end
