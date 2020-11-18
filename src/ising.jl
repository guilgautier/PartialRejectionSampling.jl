struct Ising
    "Grid dimension"
    dims::AbstractVector{Int}
    "Interaction graph"
    g::LG.SimpleGraph{Int}
    "Magnetization"
    h::Union{Real, AbstractVector{Real}}
    "Interaction"
    J::Real
end

function Ising(
    dims::AbstractVector{Int},
    periodic::Bool,
    h::Union{Real, AbstractVector{Real}},
    J::Real
)
    if isa(h, AbstractVector)
        @assert length(h) != prod(dims) "length(h) != prod(dims)"
    end
    g = LG.grid(dims; periodic=periodic)
    return Ising(dims, g, h, J)
end

function Ising(
    dims::AbstractVector{Int},
    periodic::Bool,
    J::Real
)
    return Ising(dims, periodic, 0.0, J)
end

function energy(
    ising::Ising,
    state::AbstractVector{Int}
)::Real
    @assert all(x -> x ∈ [-1, 1], state) "state values not all {-1, 1}"

    E = 0.0

    h = ising.h
    if isa(h, Real)
        E = - h * sum(state)
    else  # isa(h, AbstractVector)
        E = - h' * state
    end

    J = ising.J
    for e in LG.edges(ising.g)
        i, j = Tuple(e)
        xᵢ, xⱼ = state[i], state[j]
        E += J * xᵢ * xⱼ
    end
    return E
end

function plot(
    ising::Ising,
    state::AbstractVector{Int}
)

    pos = collect(Iterators.product(1:ising.dims[1], 1:ising.dims[2]))[:]
    locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

    col_nodes = ifelse.(
                    state .== 1,
                    Colors.colorant"gray",
                    Colors.colorant"white")

    p = GraphPlot.gplot(
            ising.g,
            locs_x,
            reverse(locs_y),
            nodefillc=col_nodes)
    #     nodelabel=LG.vertices(g),
    #     arrowlengthfrac=0.05
    #     edgestrokec=col_edges

    display(p)
end

function sample_state!(
    state::AbstractVector{T},
    ising::Ising,
    i::T;
    rng=-1
) where {T<:Int}
    rng = getRNG(rng)
    hᵢ = isa(ising.h, Real) ? ising.h : ising.h[i]
    state[i] = rand(rng) < sigmoid(hᵢ) ? 1 : -1
end

function sample_states!(
    state::AbstractVector{T},
    ising::Ising,
    indices::Set{T};
    rng=-1
) where {T<:Int}
    for i in indices
        sample_state!(state, ising, i; rng=rng)
    end
end

function sample_conditional!(
    state::AbstractVector{T},
    ising::Ising,
    i::T;
    rng=-1
) where {T<:Int}
    rng = getRNG(rng)
    hᵢ = isa(ising.h, Real) ? ising.h : ising.h[i]
    proba = sigmoid(hᵢ + ising.J * sum(state[j] for j in LG.neighbors(ising.g, i)))
    state[i] = rand(rng) < proba ? 1 : -1
end

## Partial Rejection Sampling sampler (mix of Guo Jerrum Liu and Moka Kroese)

function weighted_graph(
    ising::Ising;
    rng=-1
)::SWG.SimpleWeightedGraph{Int64,Float64}
    rng = getRNG(rng)
    wg = SWG.SimpleWeightedGraph(ising.g)
    for e in LG.edges(wg)
        i, j = Tuple(e)
        @inbounds wg.weights[i, j] = wg.weights[j, i] = rand(rng)
    end
    return wg
end

function bad_states(
    ising::Ising,
    state::AbstractVector{T},
    wg::SWG.SimpleWeightedGraph{T, U}
)::Set{T} where {T, U}

    _2absJ, sign_J = -2.0 * abs(ising.J), sign(ising.J)

    bad = Set{T}()
    for e ∈ LG.edges(wg)
        i, j = Tuple(e)
        # Break constraint:
        # log(U_ij) > - |J| (1 - sign(J) x_i x_j)
        # <=> sign(J) x_i x_j < 0 & log(U_ij) > - 2 |J|
        if ((sign_J * state[i] * state[j]) < 0.0) & (log(e.weight) > _2absJ)
            union!(bad, i, j)
        end
    end
    return bad
end

function resampling_states!(
    wg::SWG.SimpleWeightedGraph{T, U},
    ising::Ising,
    state::AbstractVector{T};
    rng=-1
)::Set{T} where {T, U}

    rng = getRNG(rng)
    sign_J = sign(ising.J)
    R = bad_states(ising, state, wg)
    ∂R, ∂R_tmp = copy(R), Set{T}()
    while !isempty(∂R)
        for i ∈ ∂R
            for j ∈ LG.neighbors(wg, i)
                # Break constraint:
                # log(U_ij) > - |J| (1 - sign(J) x_i x_j)
                # <=> sign(J) x_i x_j < 0 & log(U_ij) > - 2 |J|
                if j ∈ R
                    if (sign_J * state[i] * state[j]) < 0.0
                        # U_ij can be increased
                        @inbounds wg.weights[i, j] = wg.weights[j, i] = rand(rng)
                    end
                else # if j ∉ R
                    # x_j can be flipped to make sign(J) x_i x_j < 0
                    union!(R, j)
                    union!(∂R_tmp, j)
                    # followed by an increase of U_ij
                    @inbounds wg.weights[i, j] = wg.weights[j, i] = rand(rng)
                end
            end
        end
        ∂R, ∂R_tmp = ∂R_tmp, Set{T}()
    end
    return R
end

function prs(
    ising::Ising;
    rng=-1,
)::Tuple{AbstractVector{Int}, Int}
    rng = getRNG(rng)

    g = weighted_graph(ising; rng=rng)
    state = Vector{Int}(undef, LG.nv(g))
    res = Set{Int}(1:LG.nv(g))

    cnt = 0
    while !isempty(res)
        sample_states!(state, ising, res; rng=rng)
        res = resampling_states!(g, ising, state; rng=rng)
        cnt += 1
    end

    return state, cnt
end

## Perfect Gibbs sampler (Feng, Guo, Yin) https://arxiv.org/pdf/1907.06033.pdf

function bayes_filter(
    ising::Ising,
    state::AbstractVector{T},
    i::T,
    R::Set{T};
    rng=-1
)::Bool where {T<:Int}

    ∂i_R̄ = setdiff(LG.neighbors(ising.g, i), R)
    isempty(∂i_R̄) && return true

    rng = getRNG(rng)
    Xi_J = state[i] * ising.J
    acc_r = Xi_J * (-sign(Xi_J) * length(∂i_R̄) - sum(state[j] for j in ∂i_R̄))
    return log(rand(rng)) < acc_r
end

function gibbs_perfect(
    ising::Ising;
    rng=-1
)::Tuple{AbstractVector{Int}, Int}
    rng = getRNG(rng)

    n = LG.nv(ising.g)
    R = Set{Int}(1:n)
    state = Vector{Int}(undef, n)
    sample_states!(state, ising, R; rng=rng)

    cnt = 1
    while !isempty(R)
        i = rand(rng, R)
    if bayes_filter(ising, state, i, R; rng=-1)
            sample_conditional!(state, ising, i; rng=rng)
            delete!(R, i)
        else
            union!(R, LG.neighbors(ising.g, i))
        end
        cnt += 1
    end

    return state, cnt
end
