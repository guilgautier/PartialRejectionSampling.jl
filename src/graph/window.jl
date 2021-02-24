"""
    GraphNode{T<:Int64} <: AbstractDiscreteWindow{T}

Structure with unique field `idx` representing the index of the vertex of a graph.
"""
struct GraphNode{T<:Int64} <: AbstractDiscreteWindow{T}
    # Corner
    idx::T
end
