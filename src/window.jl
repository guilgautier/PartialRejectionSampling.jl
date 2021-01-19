abstract type AbstractWindow end

# Discrete

abstract type AbstractDiscreteWindow{T} <: AbstractWindow end

struct GraphNode{T<:Int64} <: AbstractDiscreteWindow{T}
    # Corner
    idx::T
end
function GraphNode(idx::Integer)
    return GraphNode{typeof(idx)}(idx)
end

# Spatial

abstract type AbstractSpatialWindow{T<:Float64} <: AbstractWindow end

abstract type AbstractRectangleWindow{T} <: AbstractSpatialWindow{T} end
# Πᵢ [cᵢ, cᵢ + wᵢ]
struct RectangleWindow{T} <: AbstractRectangleWindow{T}
    # Corner
    c::Vector{T}
    # Width
    w::Vector{T}
end

function RectangleWindow(c::Vector, w::Vector)
    @assert length(c) == length(w)
    return RectangleWindow{Float64}(c, w)
end

# Πᵢ [cᵢ, cᵢ + w]
struct SquareWindow{T} <: AbstractRectangleWindow{T}
    # Corner
    c::Vector{T}
    # Width
    w::T
end

function SquareWindow(c::Vector, w::Real)
    return SquareWindow{Float64}(c, w)
end

struct BallWindow{T} <: AbstractSpatialWindow{T}
    center::Vector{T}
    radius::T
end

function BallWindow(center::Vector, radius::Real)
    return BallWindow{Float64}(center, radius)
end

dimension(win::GraphNode) = 0
dimension(win::AbstractSpatialWindow) = length(win.c)

volume(win::GraphNode) = 0
volume(win::RectangleWindow) = prod(win.w)
volume(win::SquareWindow) = win.w^dimension(win)
function volume(win::BallWindow)
    d = dimension(win)
    return π^(d/2) * win.radius^d / gamma(d/2 + 1)
end

function Base.in(
    x::AbstractVector,
    win::AbstractRectangleWindow
)
    return all(win.c .<= x) && all(x .<= win.c .+ win.w)
end

function Base.in(
    x::AbstractVector,
    win::BallWindow
)
    return Distances.euclidean(x, win.center) <= win.radius
end

"""
Sample uniformly in `RectangleWindow`
"""
function Base.rand(
    win::AbstractRectangleWindow{T};
    rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    return win.c .+ win.w .* rand(rng, T, dimension(win))
end

"""
Sample uniformly in `BallWindow`
"""
function Base.rand(
    win::BallWindow{T};
    rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    d = dimension(win)
    x = zeros(T, d+2)
    randn!(rng, x)
    normalize!(x)
    return win.radius .* x[1:d] .+ win.center
end
