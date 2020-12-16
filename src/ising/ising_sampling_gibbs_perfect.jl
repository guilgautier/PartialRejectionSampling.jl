## Perfect Gibbs sampler (Feng, Guo, Yin) https://arxiv.org/pdf/1907.06033.pdf

function generate_sample_gibbs_perfect(
    ising::Ising{T};
    rng=-1
)::Vector{T} where {T<:Integer}
    rng = getRNG(rng)

    n = LG.nv(ising.g)
    R = Set{T}(1:n)
    state = Vector{T}(undef, n)
    generate_sample!(state, R, ising; rng=rng)

    while !isempty(R)
        i = rand(rng, R)
        if bayes_filter(ising, state, i, R; rng=-1)
            generate_sample_conditional!(state, i, ising; rng=rng)
            delete!(R, i)
        else
            union!(R, LG.neighbors(ising.g, i))
        end
    end

    return state
end

function generate_sample!(
        state::AbstractVector{T},
        indices,
        ising::Ising{T};
        rng=-1
) where {T<:Integer}
    rng = getRNG(rng)
    for i in indices
        hᵢ = ising.h isa Number ? ising.h : ising.h[i]
        state[i] = rand(rng) < sigmoid(hᵢ) ? 1 : -1
    end
end

function bayes_filter(
    ising::Ising{T},
    state::AbstractVector{T},
    i::T,
    R::Set{T};
    rng=-1
)::Bool where {T<:Integer}

    ∂i_R̄ = setdiff(LG.neighbors(ising.g, i), R)
    isempty(∂i_R̄) && return true

    rng = getRNG(rng)
    Xi_J = state[i] * ising.J
    acc_ratio = Xi_J * (-sign(Xi_J) * length(∂i_R̄) - sum(state[j] for j in ∂i_R̄))
    return log(rand(rng)) < acc_ratio
end

function generate_sample_conditional!(
    state::AbstractVector{T},
    i::T,
    ising::Ising{T};
    rng=-1
) where {T<:Integer}
    rng = getRNG(rng)
    hᵢ = ising.h isa Number ? ising.h : ising.h[i]
    proba = sigmoid(hᵢ + ising.J * sum(state[j] for j in LG.neighbors(ising.g, i)))
    state[i] = rand(rng) < proba ? 1 : -1
end
