struct HardCoreGraph{T<:Integer} <: AbstractGraphPointProcess{T}
    "Graph"
    g::LG.SimpleGraph{T}
    β::Float64
end

function HardCoreGraph(
        g::LG.SimpleGraph{T},
        β::Real
) where {T<:Integer}
    @assert β >= 0
    return HardCoreGraph{T}(g, float(β))
end

function generate_sample_prs(
        hcg::HardCoreGraph{T};
        rng=-1
)::Vector{T} where {T}

    proba = hcg.β / (one(hcg.β) + hcg.β)

    rng = getRNG(rng)
    adj = LG.adjacency_matrix(hcg.g)
    occupied = randsubseq(rng, LG.vertices(hcg.g), proba)
    while true
        # Check if occupied vertices form an independent set
        sub_graph = LG.SimpleGraph(adj[occupied, occupied])
        LG.ne(sub_graph) == 0 && break

        independent, resample = T[], Set{T}()
        for cc in LG.connected_components(sub_graph)
            if length(cc) == 1  # Identify current independent vertices
                append!(independent, occupied[cc])
            else  # Construct the resampling set of vertices
                union!(resample, occupied[cc])
                for v in cc
                    union!(resample, LG.neighbors(hcg.g, occupied[v]))
                end
            end
        end
        randsubseq!(rng, occupied, collect(resample), proba)
        append!(occupied, independent)
    end
    return occupied
end
