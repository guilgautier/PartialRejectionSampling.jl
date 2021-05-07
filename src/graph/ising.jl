@doc raw"""
    Ising{T<:Int} <: AbstractGraphPointProcess{T}

The [Ising model](https://en.wikipedia.org/wiki/Ising_model) characterizes a point process defined on the vertices of a graph ``(V, E)`` with joint density proportional to

```math
    \mathbb{P}\!\left[ \mathcal{X}=X \right]
    \propto
    \prod_{i \in V}
        \exp(h_i x_i)
    \prod_{\{i, j\} \in E}
        \exp(J x_i x_j)
```

where ``(h_i)_{i\in V}`` are called magnetization paremeters and ``J`` the interaction coefficient (``J \gtreqless 0`` characterizes ferro/antiferro magnetic interactions).

# Example

A realization from a ``5\times 5`` grid graph, with ``h = 0.2, J = 0.1``.

![assets/ising_grid_prs.png](assets/ising_grid_prs.png)
"""
struct Ising{T<:Int} <: AbstractGraphPointProcess{T}
    graph::LG.SimpleGraph{T}
    "Interaction"
    J::Float64
    "Magnetization"
    h::Union{Float64,Vector{Float64}}
end

function Base.show(io::IO, pp::Ising{T}) where {T}
    print(io, "Ising{$T}\n- graph = $(pp.graph)\n- J = $(pp.J) (interaction)\n- h = $(pp.h) (magnetization)")
end

"""
    Ising(
        dims::Vector{T},
        J::Real,
        h::Union{Real,Vector{Real}}=0;
        periodic::Bool=true
    ) where {T<:Int}

Construct [`PRS.Ising`](@ref) on a grid graph with dimension `dims`, interaction parameter `J` and magnetization `h`.

Periodic boundary conditions on the grid graph are set according to `periodic`.

```jldoctest; output = true
using PartialRejectionSampling

dims = [5, 5]
J, h = 0.01, 0
PRS.Ising(dims, J, h; periodic=true)

# output

Ising{Int64}
- graph = {25, 50} undirected simple Int64 graph
- J = 0.01 (interaction)
- h = 0.0 (magnetization)
```
"""
function Ising(
    dims::Vector{T},
    J::Real,
    h::Union{Real,Vector{Real}}=0;
    periodic::Bool=true
) where {T<:Int}
    g = LG.grid(dims; periodic=periodic)
    return Ising(g, J, float.(h))
end

"""
    Ising(
        graph::LG.SimpleGraph{T},
        J::Real,
        h::Union{Real,Vector{Real}}=0
    ) where {T<:Int}

Construct [`PRS.Ising`](@ref) on `graph` with interaction parameter `J` and magnetization `h`.

```jldoctest; output = true
using PartialRejectionSampling
using LightGraphs; const LG = LightGraphs

graph = LG.grid([5, 5])
J, h = 0.01, 0
PRS.Ising(graph, J, h)

# output

Ising{Int64}
- graph = {25, 40} undirected simple Int64 graph
- J = 0.01 (interaction)
- h = 0.0 (magnetization)
```
"""
function Ising(
    graph::LG.SimpleGraph{T},
    J::Real,
    h::Union{Real,Vector{Real}}=0
) where {T<:Int}
    h isa Vector && @assert length(h) == LG.nv(graph)
    return Ising{T}(graph, J, float.(h))
end

## Sampling

@doc raw"""
    generate_sample_prs(
        [rng::Random.AbstractRNG,]
        pp::Ising
    )

Generate an exact sample form `pp` using Partial Rejection Sampling (PRS), see Section 4.2 of [GuJeLi19](@cite).

Default sampler is [`PRS.generate_sample_grid_prs`](@ref).

**See also**

- [FeViYi19](@cite),
- [FeGuYi19](@cite) [`PRS.generate_sample_gibbs_perfect`](@ref).
"""
generate_sample_prs(rng::Random.AbstractRNG, pp::Ising) = generate_sample_grid_prs(rng, pp)
# generate_sample_prs(pp::Ising) = generate_sample_prs(Random.default_rng(), pp)

@doc raw"""
    generate_sample(
        [rng::Random.AbstractRNG,]
        pp::Ising{T},
        idx::T
    )::T where {T<:Int}

Generate an exact sample from the marginal distribution of state ``i=`` `idx` of `pp`.

More specifically,

```math
    x_i
        \sim
        \operatorname{Bernoulli}_{-1, 1}
            (\sigma(h_i)),
```

where ``\sigma`` denotes the [sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function).
"""
function generate_sample(
    rng::Random.AbstractRNG,
    pp::Ising{T},
    idx::T
)::T where {T<:Int}
    hᵢ = pp.h isa Real ? pp.h : pp.h[idx]
    return rand(rng) < sigmoid(hᵢ) ? one(T) : -one(T)
end

# function generate_sample(
#     pp::Ising{T},
#     idx::T
# )::T where {T<:Int}
#     return generate_sample(Random.default_rng(), pp, idx)
# end

@doc raw"""
    generate_sample!(
        rng::Random.AbstractRNG,
        state::AbstractVector{T},
        indices,
        ising::Ising{T}
    ) where {T<:Int}

Generate an exact sample from the marginal distribution of each states of the [`PRS.Ising`](@ref) model `ising` at the prescribed `indices`.
"""
function generate_sample!(
    rng::Random.AbstractRNG,
    state::AbstractVector{T},
    indices,
    ising::Ising{T}
) where {T<:Int}
    for i in indices
        state[i] = generate_sample(rng, ising, i)
    end
end

@doc raw"""
    generate_sample_conditional!(
        rng::Random.AbstractRNG,
        state::AbstractVector{T},
        i::T,
        ising::Ising{T}
    ) where {T<:Int}

Generate an exact sample from the conditional distribution of the state ``x_i`` given its neighboring states in `ising.graph`.
More specifically,
```math
    x_i \mid x_{N(i)}
        \sim
        \operatorname{Bernoulli}_{-1, 1}
            (\sigma(h_i + J \sum_{j \in N(i)} x_j)),

where ``\sigma`` denotes the [sigmoid function](https://en.wikipedia.org/wiki/Sigmoid_function).
```
"""
function generate_sample_conditional!(
    rng::Random.AbstractRNG,
    state::AbstractVector{T},
    i::T,
    ising::Ising{T}
) where {T<:Int}
    hᵢ = ising.h isa Number ? ising.h : ising.h[i]
    proba = sigmoid(hᵢ + ising.J * sum(state[j] for j in LG.neighbors(ising.graph, i)))
    state[i] = rand(rng) < proba ? 1 : -1
end

include("ising_sampling_gibbs_perfect.jl")
include("ising_sampling_grid_prs.jl")
