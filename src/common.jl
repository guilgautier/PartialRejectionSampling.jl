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
    AbstractSpatialPointProcess{T<:Vector{Float64}} <: AbstractPointProcess{T}

Abstract type encoding point processes defined on ``\mathbb{R}^d``.

Concrete instances must have a `window` field of type [`PRS.AbstractSpatialWindow`](@ref)
"""
abstract type AbstractSpatialPointProcess{T<:Vector{Float64}} <: AbstractPointProcess{T} end

"""
    window(pp::AbstractSpatialPointProcess) = pp.window
"""
window(pp::AbstractSpatialPointProcess) = pp.window

"""
    dimension(pp::AbstractSpatialPointProcess) = dimension(window(pp))

**See also**
    - [`window`](@ref)
"""
dimension(pp::AbstractSpatialPointProcess) = dimension(window(pp))

abstract type AbstractGraphPointProcess{T} <: AbstractPointProcess{T} end

# Default methods

"""
Default exact sampling method
"""
function generate_sample end

"""
Default exact sampling method using Partial Rejection Sampling (PRS, [GuJeLi19](@cite))
"""
function generate_sample_prs end
