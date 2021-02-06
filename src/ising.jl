@doc raw"""
    Ising{T<:Int} <: AbstractGraphPointProcess{T}

[Ising](https://en.wikipedia.org/wiki/Ising_model) is a point process defined on the vertices of `graph` ``=(V, E)`` with joint density proportional to

```math
    \mathbb{P}\!\left[ \mathcal{X}=X \right]
    \prod_{i \in V}
        \exp(h_i x_i)
    \prod_{\{i, j\} \subseteq E}
        \exp(J x_i x_j)
```

where ``(h_i)_{V}`` are called magnetization paremeters and ``J`` the interaction coefficient (``J \pm 0`` characterizes ferro/antiferro magnetic interactions).
"""
struct Ising{T<:Int} <: AbstractGraphPointProcess{T}
    graph::LG.SimpleGraph{T}
    "Magnetization"
    h::Union{Float64,Vector{Float64}}
    "Interaction"
    J::Float64
end

"""
    Ising(
        dims::Vector{T},
        periodic::Bool,
        h::Real,
        J::Real,
    ) where {T<:Int}

Construct [`Ising`](@ref) on a grid graph with dimension `dims` with periodic boundary conditions according to `periodic`, interaction parameter `J` and constant magnetization `h`
"""
function Ising(
    dims::Vector{T},
    periodic::Bool,
    h::Real,
    J::Real
) where {T<:Int}
    g = LG.grid(dims; periodic=periodic)
    return Ising(g, float(h), J)
end

"""
    Ising(
        dims::Vector{T},
        periodic::Bool,
        h::Vector{Real},
        J::Real
    ) where {T<:Int}

Construct [`Ising`](@ref) on a grid graph with dimension `dims` with periodic boundary conditions according to `periodic`, interaction parameter `J` and magnetization vector `h`
"""
function Ising(
    dims::Vector{T},
    periodic::Bool,
    h::Vector{Real},
    J::Real
) where {T<:Int}
    @assert length(h) == prod(dims)
    g = LG.grid(dims; periodic=periodic)
    return Ising(g, float.(h), J)
end

"""
    Ising(
        graph::LG.SimpleGraph{T},
        h::Real,
        J::Real
    ) where {T<:Int}

Construct [`Ising`](@ref) on `graph` with interaction parameter `J` and constant magnetization `h`
"""
function Ising(
    graph::LG.SimpleGraph{T},
    h::Real,
    J::Real
) where {T<:Int}
    return Ising{T}(graph, float(h), J)
end

"""
    Ising(
        graph::LG.SimpleGraph{T},
        h::Vector{Real},
        J::Real
    ) where {T<:Int}

Construct [`Ising`](@ref) on `graph` with interaction parameter `J` and magnetization vector `h`
"""
function Ising(
    graph::LG.SimpleGraph{T},
    h::Vector{Real},
    J::Real
) where {T<:Int}
    return Ising{T}(graph, float.(h), J)
end

## Sampling

include("ising/ising_sampling_gibbs_perfect.jl")
include("ising/ising_sampling_grid_prs.jl")
