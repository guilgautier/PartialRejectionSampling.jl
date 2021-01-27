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

abstract type AbstractRectangleWindow{T<:Float64} <: AbstractSpatialWindow{T} end
# Πᵢ [cᵢ, cᵢ + wᵢ]
struct RectangleWindow{T<:Float64} <: AbstractRectangleWindow{T}
    # Corner
    c::Vector{T}
    # Width
    w::Vector{T}
end

function RectangleWindow(c::AbstractVector, w::Vector)
    @assert length(c) == length(w)
    @assert all(w .> 0)
    return RectangleWindow{Float64}(c, w)
end

# Πᵢ [cᵢ, cᵢ + w]
struct SquareWindow{T<:Float64} <: AbstractRectangleWindow{T}
    # Corner
    c::Vector{T}
    # Width
    w::T
end

function SquareWindow(c::AbstractVector, w::Real=1.0)
    @assert w > 0
    return SquareWindow{Float64}(c, w)
end

function rectangle_square_window(c::AbstractVector, w::Union{Real, AbstractVector})
    return w isa Real ? SquareWindow(c, w) : RectangleWindow(c, w)
end

struct BallWindow{T<:Float64} <: AbstractSpatialWindow{T}
    # Center
    c::Vector{T}
    # Radius
    r::T
end

function BallWindow(c::AbstractVector, r::Real)
    @assert r > 0
    return BallWindow{Float64}(c, r)
end

dimension(win::GraphNode) = 0
dimension(win::AbstractSpatialWindow) = length(win.c)

volume(win::GraphNode) = 0
volume(win::RectangleWindow) = prod(win.w)
volume(win::SquareWindow) = win.w^dimension(win)
function volume(win::BallWindow)
    d = dimension(win)
    return π^(d/2) * win.r^d / SF.gamma(d/2 + 1)
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
    return Distances.euclidean(x, win.c) <= win.r
end

"""
Sample uniformly in `AbstractRectangleWindow`
"""
function Base.rand(
    win::AbstractRectangleWindow{T};
    rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    return win.w .* rand(rng, dimension(win)) .+ win.c
end

"""
Sample n points uniformly in `AbstractRectangleWindow`
"""
function Base.rand(
    win::AbstractRectangleWindow{T},
    n::Int;
    rng=-1
)::Matrix{T} where {T}
    rng = getRNG(rng)
    x = Matrix{T}(undef, dimension(win), n)
    rand!(rng, x)
    return win.w .* x .+ win.c
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
    x = Vector{T}(undef, d+2)
    randn!(rng, x)
    normalize!(x)
    return win.r .* x[1:d] .+ win.c
end

"""
Sample n points uniformly in `BallWindow`
"""
function Base.rand(
    win::BallWindow{T},
    n::Int;
    rng=-1
)::Matrix{T} where {T}
    rng = getRNG(rng)
    d = dimension(win)
    X = Matrix{T}(undef, d+2, n)
    Random.randn!(rng, X)
    normalize_columns!(X, 2)
    return win.r .* X[1:d, :] .+ win.c
end
