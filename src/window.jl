abstract type AbstractWindow end

abstract type AbstractDiscreteWindow{T} <: AbstractWindow end

struct GraphNode{T<:Int64} <: AbstractDiscreteWindow{T}
    # Corner
    idx::T
end

abstract type AbstractSpatialWindow{T<:Float64} <: AbstractWindow end
# Πᵢ [cᵢ, cᵢ + wᵢ]
struct RectangleWindow{T} <: AbstractSpatialWindow{T}
    # Corner
    c::Vector{T}
    # Width
    w::Vector{T}
end

# Πᵢ [cᵢ, cᵢ + w]
struct SquareWindow{T} <: AbstractSpatialWindow{T}
    # Corner
    c::Vector{T}
    # Width
    w::T
end

function GraphNode(idx::Integer)
    return GraphNode{typeof(idx)}(idx)
end

function RectangleWindow(c::Vector, w::Vector)
    @assert length(c) == length(w)
    return RectangleWindow{Float64}(float.(c), float.(w))
end

function SquareWindow(c::Vector, w::Real)
    return SquareWindow{Float64}(float.(c), float(w))
end

function spatial_window(c::Vector, w::Union{Real,Vector})
    return allequal(w) ? SquareWindow(c, w[1]) : RectangleWindow(c, w)
end

dimension(win::GraphNode) = 0
dimension(win::AbstractSpatialWindow) = length(win.c)

volume(win::GraphNode) = 0
volume(win::RectangleWindow) = prod(win.w)
volume(win::SquareWindow) = win.w^dimension(win)

function Base.rand(
    win::AbstractSpatialWindow{T};
    rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    return win.c .+ win.w .* rand(rng, T, dimension(win))
end
