@doc raw"""
    HomogeneousPoissonPointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}

Homegeneous [Poisson point process](https://en.wikipedia.org/wiki/Poisson_point_process) with intensity ``\beta > 0``, denoted ``\operatorname{Poisson}(\beta)``.

``\operatorname{Poisson}(\beta)`` has density (w.r.t. the homogenous Poisson point process with unit intensity ``\operatorname{Poisson}(1)``) proportional to

```math
    \prod_{x \in X}
        \beta.
```
"""
struct HomogeneousPoissonPointProcess{T<:Vector{Float64}} <: AbstractSpatialPointProcess{T}
    β::Float64
    window::AbstractSpatialWindow{Float64}
end

function Base.show(io::IO, pp::HomogeneousPoissonPointProcess{T}) where {T}
    print(io, "HomogeneousPoissonPointProcess{$T}\n- β = $(pp.β)\n- window = $(pp.window)")
end

@doc raw"""
    HomogeneousPoissonPointProcess(
        β::Real,
        window::Union{Nothin,AbstractSpatialWindow}=nothing
    )

Construct a [`PRS.HomogeneousPoissonPointProcess`](@ref) with intensity `β` restricted to `window`.

Default window (`window=nothing`) is [`PRS.SquareWindow`](@ref)`()`.

```jldoctest; output = true
using PartialRejectionSampling

β = 40
win = PRS.SquareWindow(zeros(2), 1)
PRS.HomogeneousPoissonPointProcess(β, win)

# output

HomogeneousPoissonPointProcess{Array{Float64,1}}
- β = 40.0
- window = SquareWindow [0.0, 1.0]^2
```
"""
function HomogeneousPoissonPointProcess(
    β::Real,
    window::Union{Nothing,AbstractSpatialWindow}=nothing
)
    @assert β > 0
    win = window === nothing ? SquareWindow() : window
    return HomogeneousPoissonPointProcess{Vector{Float64}}(β, win)
end

"""
    intensity(pp::HomogeneousPoissonPointProcess) = pp.β
"""
intensity(pp::HomogeneousPoissonPointProcess) = pp.β

## Sampling

"""
    generate_sample(
        pp::HomogeneousPoissonPointProcess{Vector{T}};
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
    )::Matrix{T} where {T<:Float64}

Generate an exact sample from [`PRS.HomogeneousPoissonPointProcess`](@ref) on window `win`.
Sampled points are stored as columns of the output matrix.

Default window (`win=nothing`) is `window(pp)`.
"""
function generate_sample(
        pp::HomogeneousPoissonPointProcess{Vector{T}};
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
)::Matrix{T} where {T<:Float64}
    rng = getRNG(rng)
    window_ = win === nothing ? window(pp) : win
    n = rand(rng, Distributions.Poisson(pp.β * volume(window_)))
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

Generate an exact sample from a homogenous [`PRS.HomogeneousPoissonPointProcess`](@ref) Poisson(β) on ``\bigcup_{i} B(c_i, r)`` (union of balls centered at ``c_i`` with the same radius ``r``).

If `win ≂̸ nothing`, returns the points falling in `win`.

!!! hint

    Use the independence property of the Poisson point process on disjoint subsets in order to

    - Sample from Poisson(β) on ``B(c_1, r)``,
    - Sample from Poisson(β) on ``B(c_2, r) \setminus B(c_1, r)``,
    - Sample from Poisson(β) on ``B(c_j, r) \setminus \bigcup_{i<j} B(c_i, r)``,
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
