"""
Abstract type representing a window

**See also**

- [`PRS.AbstractDiscreteWindow`](@ref)
- [`PRS.AbstractSpatialWindow`](@ref)
"""
abstract type AbstractWindow end

# Discrete

"""
Abstract type representing a window on a discrete state space

**See also**

- [`PRS.GraphNode`](@ref)
"""
abstract type AbstractDiscreteWindow{T} <: AbstractWindow end

"""
    GraphNode{T<:Int64} <: AbstractDiscreteWindow{T}

Structure with unique field `idx` representing the index of the vertex of a graph
"""
struct GraphNode{T<:Int64} <: AbstractDiscreteWindow{T}
    # Corner
    idx::T
end

"""
    GraphNode(idx::T) where {T<:Int}

Construct a [`GraphNode`](@ref)
"""
function GraphNode(idx::T) where {T<:Int}
    return GraphNode{T}(idx)
end

# Spatial

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

@doc raw"""
    AbstractRectangleWindow{T<:Float64} <: AbstractSpatialWindow{T}

Abstract type representing a [hyperrectangle](https://en.wikipedia.org/wiki/Hyperrectangle) ``\prod_i [c_i, c_i + w_i]``

**See also**

- [`PRS.RectangleWindow`](@ref)
- [`PRS.SquareWindow`](@ref)
"""
abstract type AbstractRectangleWindow{T<:Float64} <: AbstractSpatialWindow{T} end

@doc raw"""
    RectangleWindow{T<:Float64} <: AbstractRectangleWindow{T}

Structure representing a hyperrectangle ``\prod_i [c_i, c_i + w_i]``, with fields
- `c` lower left corner of the hyperrectangle
- `w` width vector of the hyperrectangle along each coordinate
"""
struct RectangleWindow{T<:Float64} <: AbstractRectangleWindow{T}
    # Corner
    c::Vector{T}
    # Width
    w::Vector{T}
end

"""
    RectangleWindow(c::AbstractVector, w::Vector)

Construct a [`PRS.RectangleWindow`](@ref)
"""
function RectangleWindow(c::AbstractVector, w::Vector)
    @assert length(c) == length(w)
    @assert all(w .> 0)
    return RectangleWindow{Float64}(c, w)
end

@doc raw"""
    SquareWindow{T<:Float64} <: AbstractRectangleWindow{T}

Structure representing a [hypercube](https://fr.wikipedia.org/wiki/Hypercube) ``\prod_i [c_i, c_i + w]``, with fields
- `c` lower left corner of the hypercube
- `w` length of the hypercube
"""
struct SquareWindow{T<:Float64} <: AbstractRectangleWindow{T}
    # Corner
    c::Vector{T}
    # Width
    w::T
end

"""
    SquareWindow(c::AbstractVector, w::Real=1.0)

Construct a [`PRS.SquareWindow`](@ref)
"""
function SquareWindow(c::AbstractVector, w::Real=1.0)
    @assert w > 0
    return SquareWindow{Float64}(c, w)
end

"""
    rectangle_square_window(c, w)

Construct a [`PRS.RectangleWindow`](@ref) or a [`PRS.SquareWindow`](@ref) depending on whether all coordinates of `w` are the same
"""
function rectangle_square_window(c, w)
    return allequal(w) ? SquareWindow(c, w[1]) : RectangleWindow(c, w)
end

@doc raw"""
    BallWindow{T<:Float64} <: AbstractSpatialWindow{T}

Structure representing a closed [ball](https://en.wikipedia.org/wiki/Ball_(mathematics)) ``B(c, r)``, with fields
- `c` center of the ball
- `r` radius of the ball
"""
struct BallWindow{T<:Float64} <: AbstractSpatialWindow{T}
    # Center
    c::Vector{T}
    # Radius
    r::T
end

"""
    BallWindow(c::AbstractVector, r::Real)

Construct a [`PRS.BallWindow`](@ref)
"""
function BallWindow(c::AbstractVector, r::Real=1.0)
    @assert r > 0
    return BallWindow{Float64}(c, r)
end

"""
    dimension(win::GraphNode) = 0
"""
dimension(win::GraphNode) = 0

"""
    dimension(win::AbstractSpatialWindow) = length(win.c)

Return the dimension of window `win`
"""
dimension(win::AbstractSpatialWindow) = length(win.c)

"""
    volume(win::GraphNode) = 0
"""
volume(win::GraphNode) = 0

@doc raw"""
    volume(win::RectangleWindow) = prod(win.w)
"""
volume(win::RectangleWindow) = prod(win.w)

@doc raw"""
    volume(win::SquareWindow) = win.w^dimension(win)
"""
volume(win::SquareWindow) = win.w^dimension(win)

@doc raw"""
    volume(win::BallWindow) =

Return the [volume of the ball](https://en.wikipedia.org/wiki/Volume_of_an_n-ball) ``B(c, r)\subseteq R^d`` as

```math
    \frac{π^{d/2} r^d}{\Gamma(d/2 + 1)}
```
"""
function volume(win::BallWindow)
    d = dimension(win)
    return π^(d/2) * win.r^d / SF.gamma(d/2 + 1)
end

@doc raw"""
    Base.in(
        x::AbstractVector,
        win::AbstractRectangleWindow
    )

Check if ``x \in \prod_{i} [c_i, c_i + w_i]``
"""
function Base.in(
    x::AbstractVector,
    win::AbstractRectangleWindow
)
    return all(win.c .<= x) && all(x .<= win.c .+ win.w)
end

@doc raw"""
    Base.in(
        x::AbstractVector,
        win::AbstractRectangleWindow
    )

Check if ``x \in B(c, r)``, i.e., ``\left\| x - c \right\| \leq r``
"""
function Base.in(
    x::AbstractVector,
    win::BallWindow
)
    return Distances.euclidean(x, win.c) <= win.r
end

"""
    Base.rand(
        win::AbstractRectangleWindow{T};
        rng=-1
    )::Vector{T} where {T}

Sample uniformly at random in window `win`
"""
function Base.rand(
    win::AbstractRectangleWindow{T};
    rng=-1
)::Vector{T} where {T}
    rng = getRNG(rng)
    return win.w .* rand(rng, dimension(win)) .+ win.c
end

"""
    Base.rand(
        win::AbstractRectangleWindow{T},
        n::Int;
        rng=-1
    )::Vector{T} where {T}

Sample `n` points uniformly at random in window `win`
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
    Base.rand(
        win::BallWindow{T};
        rng=-1
    )::Vector{T} where {T}

Sample uniformly at random in window `win`
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
    Base.rand(
        win::BallWindow{T},
        n::Int;
        rng=-1
    )::Vector{T} where {T}

Sample `n` points uniformly at random in window `win`
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
