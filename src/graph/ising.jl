@doc raw"""
    Ising{T<:Int} <: AbstractGraphPointProcess{T}

The [Ising model](https://en.wikipedia.org/wiki/Ising_model) characterizes a point process defined on the vertices of graph ``(V, E)`` with joint density proportional to

```math
    \mathbb{P}\!\left[ \mathcal{X}=X \right]
    \propto
    \prod_{i \in V}
        \exp(h_i x_i)
    \prod_{\{i, j\} \in E}
        \exp(J x_i x_j)
```

where ``(h_i)_{V}`` are called magnetization paremeters and ``J`` the interaction coefficient (``J \gtreqless 0`` characterizes ferro/antiferro magnetic interactions).

# Example

A realization from a ``5\times 5`` grid graph, with ``h = 0.2, J = 0.1``.

![assets/ising.png](assets/ising.png)
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

@doc raw"""
    generate_sample_prs(pp::Ising; rng=-1)

Generate an exact sample form `pp` using Partial Rejection Sampling (PRS), see Section 4.2 of [GuJeLi19](@cite)

Default sampler is [`PRS.generate_sample_grid_prs`](@ref)

**See also**

- [FeViYi19](@cite)
- [FeGuYi19](@cite) [`PRS.generate_sample_gibbs_perfect`](@ref)
"""
generate_sample_prs(pp::Ising; rng=-1) = generate_sample_grid_prs(pp; rng=rng)

@doc raw"""
    generate_sample(
        pp::Ising{T},
        idx::T;
        rng=-1
    )::T where {T<:Int}

Generate an exact sample from the marginal distribution of state ``i=`` `idx` of `pp`.

More specifically,

```math
    x_i
        \sim
        \operatorname{Bernoulli}_{-1, 1}
            (\sigma(h_i)),
```

where ``\sigma`` denotes the [`PRS.sigmoid`](@ref) function.
"""
function generate_sample(
    pp::Ising{T},
    idx::T;
    rng=-1
)::T where {T<:Int}
    rng = getRNG(rng)
    hᵢ = pp.h isa Real ? pp.h : pp.h[idx]
    return rand(rng) < sigmoid(hᵢ) ? one(T) : -one(T)
end

## Sampling

include("ising_sampling_gibbs_perfect.jl")
include("ising_sampling_grid_prs.jl")
