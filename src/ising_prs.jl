import Random
getRNG(seed::Integer = -1) = seed >= 0 ? Random.MersenneTwister(seed) : Random.GLOBAL_RNG
getRNG(seed::Union{Random.MersenneTwister,Random._GLOBAL_RNG}) = seed

import LightGraphs
const LG = LightGraphs

import SimpleWeightedGraphs
const SWG = SimpleWeightedGraphs

using GraphPlot, Colors

function ising_dependency_graph(
    dims::Tuple{T,T};
    periodic::Bool = false,
    rng = -1,
)::SWG.SimpleWeightedGraph{Int64,Float64} where {T}
    rng = getRNG(rng)
    g = SWG.SimpleWeightedGraph(LG.grid(dims, periodic = periodic))
    for e in LG.edges(g)
        i, j = Tuple(e)
        @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng)
    end
    return g
end

"""
Find bad
U_ij > exp(- |J| (1 - sign(J) x_i x_j))
log(U_ij) > - |J| (1 - sign(J) x_i x_j)
sign(J) x_i x_j < 0 & log(U_ij) > - 2 |J|
"""
function ising_find_bad_states(
    g::SWG.SimpleWeightedGraph{T, U},
    states::Vector{T},
    J::U;
    rng = -1,
)::Set{T} where {T, U}

    @assert LG.nv(g) == length(states)

    rng = getRNG(rng)

    _2absJ, sign_J = -2.0 * abs(J), sign(J)

    bad = Set{T}()
    for e ∈ LG.edges(g)
        i, j = Tuple(e)
        # Break constraint:
        # log(U_ij) > - |J| (1 - sign(J) x_i x_j)
        # <=> sign(J) x_i x_j < 0 & log(U_ij) > - 2 |J|
        if ((sign_J * states[i] * states[j]) < 0.0) & (log(e.weight) > _2absJ)
            union!(bad, i, j)
        end
    end
    return bad
end

"""
Find Res
"""
function ising_find_states_to_resample!(
    g::SWG.SimpleWeightedGraph{T, U},
    states::Vector{T},
    J::U;
    rng = -1
)::Set{T} where {T, U}

    rng = getRNG(rng)
    sign_J = sign(J)
    R = ising_find_bad_states(g, states, J, rng = rng)
    ∂R, ∂R_tmp = copy(R), Set{T}()
    while !isempty(∂R)
        for i ∈ ∂R
            for j ∈ LG.neighbors(g, i)
                # Break constraint:
                # log(U_ij) > - |J| (1 - sign(J) x_i x_j)
                # <=> sign(J) x_i x_j < 0 & log(U_ij) > - 2 |J|
                if j ∈ R
                    if (sign_J * states[i] * states[j]) < 0.0
                        # U_ij can be increased
                        @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng)
                    end
                else # if j ∉ R
                    # x_j can be flipped to make sign(J) x_i x_j < 0
                    union!(R, j)
                    union!(∂R_tmp, j)
                    # followed by an increase of U_ij
                    @inbounds g.weights[i, j] = g.weights[j, i] = rand(rng)
                end
            end
        end
        ∂R, ∂R_tmp = ∂R_tmp, Set{T}()
    end
    return R
end

"""
Resample
"""
function ising_sample_states!(
    states::Vector{T},
    res_ind::Set{T},
    probas::Vector{Float64};
    rng = -1,
) where {T}
    n = length(states)
    @assert n == length(probas)
    rng = getRNG(rng)
    for i in res_ind
        states[i] = rand(rng) < probas[i] ? 1 : -1
    end
end

sigmoid(x) = @. 1 / (1 + exp(-x))

function ising_prs(
    dims::Tuple{T, T},
    h::Vector{U},
    J::U;
    periodic::Bool=false,
    rng=-1,
) where {T, U}
    rng = getRNG(rng)

    g = ising_dependency_graph(dims, periodic=periodic, rng=rng)

    n = LG.nv(g)
    states = Vector{Int64}(undef, n)
    probas = sigmoid.(2.0 .* h)
    res = Set{T}(1:LG.nv(g))

    cnt = 0
    while !isempty(res)
        ising_sample_states!(states, res, probas, rng=rng)
        res = ising_find_states_to_resample!(g, states, J, rng = rng)
        cnt += 1
    end

    return g, states, cnt
end

function ising_prs(
    dims::Tuple{T, T},
    h::U,
    J::U;
    periodic::Bool = false,
    rng = -1,
) where {T, U}
    h_vec = fill(h, prod(dims))
    return ising_prs(dims, h_vec, J, periodic=periodic, rng=rng)
end

function plot_ising(g, dims, state)
    pos = collect(Iterators.product(1:dims[1], 1:dims[2]))[:]
    locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

    col_nodes = ifelse.(state .== 1, colorant"gray", colorant"white")

    p = gplot(g,
        locs_x, reverse(locs_y),
        nodefillc=col_nodes,
    #     nodelabel=LG.vertices(g),
    #     arrowlengthfrac=0.05
    #     edgestrokec=col_edges
    )
    display(p)
end


dims = (14, 14) # if > (14, 14) the display becomes all black, don't know why !
H, J = 0.0, -0.02 # Use Float

periodic = false
seed = -1
g, config, cnt = ising_prs(dims, H, J; periodic=periodic, rng=seed)

plot_ising(g, dims, config)

println("Number of resampling steps")
println(cnt)
