struct HomogeneousPoissonPointProcess{T} <: AbstractSpatialPointProcess{T}
    β::Float64
    window::AbstractSpatialWindow{Float64}
end

function HomogeneousPoissonPointProcess(β::Real, window::AbstractSpatialWindow)
    @assert β > 0
    return HomogeneousPoissonPointProcess{Vector{Float64}}(β, window)
end

## Sampling

"""
Sample from a homogenous Poisson point process `hp` on window `win` for which the `volume` function is defined.
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

"""
Sample from a homogenous Poisson point process Poisson(β) on ⋃ᵢ B(cᵢ, r) (union of balls centered at cᵢ with the same radius r)

Use the independence property of the Poisson point process on disjoint subsets in order to

  - Sample from Poisson(β) on B(c₁, r)
  - Sample from Poisson(β) on B(c₂, r) ∖ B(c₁, r)
  - Sample from Poisson(β) on B(cⱼ, r) ∖ ⋃_i<j B(cᵢ, r)
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
