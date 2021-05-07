"""
    AbstractPointProcess{T}

Abstract type encoding the notion of [point process](https://en.wikipedia.org/wiki/Point_process)
The type `T` corresponds to the elements' type of the point process (vectors, graph edges, etc.)

**See also**
- [`PRS.AbstractSpatialPointProcess`](@ref)
- [`PRS.AbstractGraphPointProcess`](@ref)
"""
abstract type AbstractPointProcess{T} end
Base.eltype(pp::AbstractPointProcess{T}) where {T} = T


@doc raw"""
    AbstractGraphPointProcess{T<:Vector{Float64}} <: AbstractPointProcess{T}

Abstract type encoding point processes defined on graphs
"""
abstract type AbstractGraphPointProcess{T} <: AbstractPointProcess{T} end

@doc raw"""
    AbstractSpatialPointProcess{T<:Vector{Float64}} <: AbstractPointProcess{T}

Abstract type encoding point processes defined on ``\mathbb{R}^d``.

Concrete instances must have a `window` field of type [`PRS.AbstractSpatialWindow`](@ref)
"""
abstract type AbstractSpatialPointProcess{T<:Vector{Float64}} <: AbstractPointProcess{T} end

"""
    window(pp::AbstractSpatialPointProcess)::AbstractSpatialWindow = pp.window
"""
function window(pp::AbstractSpatialPointProcess)::AbstractSpatialWindow
    return pp.window
end

"""
    dimension(pp::AbstractSpatialPointProcess) = dimension(window(pp))
"""
dimension(pp::AbstractSpatialPointProcess) = dimension(window(pp))

"""
Abstract type representing a window

**See also**

- [`PRS.AbstractDiscreteWindow`](@ref)
- [`PRS.AbstractSpatialWindow`](@ref)
"""
abstract type AbstractWindow end

@doc raw"""
    AbstractSpatialWindow{T<:Float64} <: AbstractWindow{T}

Abstract type representing a spatial window ``\subseteq \mathbb{R}^d``

**See also**

- [`PRS.AbstractRectangleWindow`](@ref)
    - [`PRS.RectangleWindow`](@ref)
    - [`PRS.SquareWindow`](@ref)
- [`PRS.BallWindow`](@ref)
"""
abstract type AbstractSpatialWindow{T<:Float64} <: AbstractWindow end

"""
Abstract type representing a window on a discrete state space

**See also**

- [`PRS.GraphNode`](@ref)
"""
abstract type AbstractDiscreteWindow{T} <: AbstractWindow end

function generate_sample(pp::AbstractPointProcess, args...; kwargs...)
    return generate_sample(Random.default_rng(), pp, args...; kwargs...)
end

function generate_sample_prs(pp::AbstractPointProcess, args...; kwargs...)
    return generate_sample_prs(Random.default_rng(), pp, args...; kwargs...)
end
