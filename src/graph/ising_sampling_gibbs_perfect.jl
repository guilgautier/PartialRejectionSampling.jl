
"""
    generate_sample_gibbs_perfect(
        [rng::Random.AbstractRNG,]
        ising::Ising{T}
    )::Vector{T} where {T<:Int}

Generate an exact realization of the [`PRS.Ising`](@ref) model using a tailored implementation of the perfect Gibbs sampler of [FeGuYi19](@cite).
"""
function generate_sample_gibbs_perfect(
    rng::Random.AbstractRNG,
    ising::Ising{T}
)::Vector{T} where {T<:Int}

    n = LG.nv(ising.graph)
    R = Set{T}(1:n)
    state = Vector{T}(undef, n)
    generate_sample!(rng, state, R, ising)

    while !isempty(R)
        i = rand(rng, R)
        if bayes_filter(rng, ising, state, i, R)
            generate_sample_conditional!(rng, state, i, ising)
            delete!(R, i)
        else
            union!(R, LG.neighbors(ising.graph, i))
        end
    end

    return state
end

function generate_sample_gibbs_perfect(ising::Ising)
    return generate_sample_gibbs_perfect(Random.default_rng(), ising)
end

"""
    bayes_filter(
        rng::Random.AbstractRNG,
        ising::Ising{T},
        state::AbstractVector{T},
        i::T,
        R::Set{T}
    )::Bool

This function is used as a subroutine of [`PRS.generate_sample_gibbs_perfect`](@ref).
"""
function bayes_filter(
    rng::Random.AbstractRNG,
    ising::Ising{T},
    state::AbstractVector{T},
    i::T,
    R::Set{T}
)::Bool where {T<:Int}
    ∂i_R̄ = setdiff(LG.neighbors(ising.graph, i), R)
    isempty(∂i_R̄) && return true

    Xi_J = state[i] * ising.J
    log_acc_ratio = Xi_J * (-sign(Xi_J) * length(∂i_R̄) - sum(state[j] for j in ∂i_R̄))
    return log(rand(rng)) < log_acc_ratio
end
