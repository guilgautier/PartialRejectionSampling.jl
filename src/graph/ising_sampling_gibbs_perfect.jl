
"""
    generate_sample_gibbs_perfect(
        ising::Ising{T};
        rng=-1
    )::Vector{T} where {T<:Int}

Generate an exact realization of the [`PRS.Ising`](@ref) model using a tailored implementation of the perfect Gibbs sampler of [FeGuYi19](@cite).
"""
function generate_sample_gibbs_perfect(
    ising::Ising{T};
    rng=-1
)::Vector{T} where {T<:Int}
    rng = getRNG(rng)

    n = LG.nv(ising.graph)
    R = Set{T}(1:n)
    state = Vector{T}(undef, n)
    generate_sample!(state, R, ising; rng=rng)

    while !isempty(R)
        i = rand(rng, R)
        if bayes_filter(ising, state, i, R; rng=rng)
            generate_sample_conditional!(state, i, ising; rng=rng)
            delete!(R, i)
        else
            union!(R, LG.neighbors(ising.graph, i))
        end
    end

    return state
end

@doc raw"""
    generate_sample!(
        state::AbstractVector{T},
        indices,
        ising::Ising{T};
        rng=-1
    ) where {T<:Int}

Generate an exact sample from the marginal distribution of each states of the [`Ising`](@ref) model `ising` at the prescribed `indices`.
"""
function generate_sample!(
    state::AbstractVector{T},
    indices,
    ising::Ising{T};
    rng=-1
) where {T<:Int}
    for i in indices
        state[i] = generate_sample(ising, i; rng=rng)
    end
end

"""
    bayes_filter(
        ising::Ising{T},
        state::AbstractVector{T},
        i::T,
        R::Set{T};
        rng=-1
    )::Bool

This function is used as a subroutine of [`PRS.generate_sample_gibbs_perfect`](@ref).
"""
function bayes_filter(
    ising::Ising{T},
    state::AbstractVector{T},
    i::T,
    R::Set{T};
    rng=-1
)::Bool where {T<:Int}

    ∂i_R̄ = setdiff(LG.neighbors(ising.graph, i), R)
    isempty(∂i_R̄) && return true

    rng = getRNG(rng)
    Xi_J = state[i] * ising.J
    acc_ratio = Xi_J * (-sign(Xi_J) * length(∂i_R̄) - sum(state[j] for j in ∂i_R̄))
    return log(rand(rng)) < acc_ratio
end

@doc raw"""
    generate_sample_conditional!(
        state::AbstractVector{T},
        i::T,
        ising::Ising{T};
        rng=-1
    ) where {T<:Int}

Generate an exact sample from the conditional distribution of the state ``x_i`` given its neighboring states in `ising.graph`.
More specifically,
```math
    x_i \mid x_{N(i)}
        \sim
        \operatorname{Bernoulli}_{-1, 1}
            (\sigma(h_i + J \sum_{j \in N(i)} x_j).
```
"""
function generate_sample_conditional!(
    state::AbstractVector{T},
    i::T,
    ising::Ising{T};
    rng=-1
) where {T<:Int}
    rng = getRNG(rng)
    hᵢ = ising.h isa Number ? ising.h : ising.h[i]
    proba = sigmoid(hᵢ + ising.J * sum(state[j] for j in LG.neighbors(ising.graph, i)))
    state[i] = rand(rng) < proba ? 1 : -1
end
