using Plots
using GraphPlot, Colors

using LazySets
const LS = LazySets

# Pedagogy
function pedagogy_generate_sample_prs(
        hc::HardCorePointProcess{T};
        win::Union{Nothing,AbstractWindow}=nothing,
        path="",
        rng=-1
)::Vector{T} where {T}

    rng = getRNG(rng)
    window_ = win === nothing ? window(hc) : win

    n = rand(rng, Distributions.Poisson(hc.β * volume(window_)))
    points = Matrix{Float64}(undef, dimension(hc), n)
    for x in eachcol(points)
        x .= rand(window_; rng=rng)
    end

    i = 0
    while true

        p = pedagogy_plot(points, true, 0.0, "white", hc.window)
        Plots.savefig(p, "$(path)/hard_core_$(i)_0.pdf")

        p = pedagogy_plot(points, true, hc.r/2, "white", hc.window)
        Plots.savefig(p, "$(path)/hard_core_$(i)_1.pdf")

        bad = vec(any(pairwise_distances(points) .< hc.r, dims=2))
        !any(bad) && break

        pedagogy_plot!(p, points[:, bad], false, hc.r/2, "red", hc.window)
        Plots.savefig(p, "$(path)/hard_core_$(i)_1bad_region.pdf")

        pedagogy_plot!(p, points[:, bad], false, hc.r, "orange", hc.window)
        Plots.savefig(p, "$(path)/hard_core_$(i)_2resample_region.pdf")

        p = pedagogy_plot(points[:, .!bad], true, 0.0, "white", hc.window)
        pedagogy_plot!(p, points[:, bad], false, hc.r, "orange", hc.window)
        Plots.savefig(p, "$(path)/hard_core_$(i)_3resample_region.pdf")

        resampled = generate_sample_poisson_union_balls(hc.β, points[:, bad], hc.r; win=window_, rng=rng)
        pedagogy_plot!(p, resampled, true, 0, "blue", hc.window)
        Plots.savefig(p, "$(path)/hard_core_$(i)_4resampled_points.pdf")

        points = hcat(points[:, .!bad], resampled)
        i += 1

    end
    return [points[:, i] for i in 1:size(points, 2)]
end

function pedagogy_plot!(
        p,
        points,
        show_center=true,
        radius=0.0,
        color="white",
        window=SquareWindow(zeros(2), 1.0)
)
    for x in (points isa Matrix ? eachcol(points) : points)
        if show_center
            Plots.scatter!(p, [x[1]], [x[2]], markersize=2, color=color)
        end
        if radius > 0
            Plots.plot!(p, LS.Ball2(vec(x), radius), color=color)
        end
    end

    Plots.xlims!(window.c[1], window.c[1] + window.w)
    Plots.ylims!(window.c[2], window.c[2] + window.w)

    return p
end

function pedagogy_plot(
        points,
        show_center=true,
        radius=0.0,
        color="white",
        window=SquareWindow(zeros(2), 1.0)
)

    p = Plots.plot([0], [0],
            label="", legend=false,
            color="white",
            linewidth=0.0,
            aspect_ratio=:equal,
            grid=:none,
            title="")

    pedagogy_plot!(p, points, show_center, radius, color, window)

    return p
end

function plot(
        pp::AbstractSpatialPointProcess,
        points;
        title=""
)
    p = Plots.plot([0], [0],
            label="", legend=false,
            color="white",
            linewidth=0.0,
            aspect_ratio=:equal,
            grid=:none,
            title=title)

    θ = collect(range(0, 2π, length=15))
    rad = pp.r / 2  # radius = interaction range / 2
    circ_x, circ_y = rad .* cos.(θ), rad .* sin.(θ)

    for x in points
        Plots.plot!(x[1] .+ circ_x,
                    x[2] .+ circ_y,
                    color="black",
                    linewidth=0.5)
    end

    win = pp.window
    Plots.xlims!(win.c[1], win.c[1] + win.w[1])
    Plots.ylims!(win.c[2], win.c[2] + (win.w isa Number ? win.w[1] : win.w[2]))

    return p
end

function plot(
        ising::Ising,
        state,
        width::Int=0,
        height::Int=0
)
    if width == 0 || height == 0
        p = GraphPlot.gplot(ising.graph,
                nodelabel=LG.vertices(ising.graph),
                nodefillc=col_nodes
                )
        return p
    else
        pos = collect(Iterators.product(1:height, 1:width))[:]
        locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

        col_nodes = ifelse.(
                    state .== 1,
                    Colors.colorant"gray",
                    Colors.colorant"white")

        p = GraphPlot.gplot(
                ising.graph,
                locs_x,
                reverse(locs_y),
                nodefillc=col_nodes)
        return p
    end
end

function plot(
    hcg::HardCoreGraph,
    state,
    width::Int=0,
    height::Int=0
)
    col_nodes = [Colors.colorant"turquoise" for _ in 1:LG.nv(hcg.graph)]
    col_nodes[state] .= Colors.colorant"red"

    if width == 0 || height == 0
        p = GraphPlot.gplot(hcg.graph,
                nodelabel=LG.vertices(hcg.graph),
                nodefillc=col_nodes
                )
        return p
    else
        pos = collect(Iterators.product(1:height, 1:width))[:]
        locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

        p = GraphPlot.gplot(hcg.graph,
                locs_x,
                reverse(locs_y),
                nodelabel=LG.vertices(hcg.graph),
                nodefillc=col_nodes
                )
        return p
    end
end

function plot(
        rsf::RootedSpanningForest,
        sample,
        width::Int=0,
        height::Int=0
)
    col_nodes, col_edges = color_cycles(sample)
    col_nodes[collect(rsf.roots)] .= Colors.colorant"orange"

    if width == 0 || height == 0
        p = GraphPlot.gplot(sample,
            nodelabel=LG.vertices(sample),
            nodefillc=col_nodes,
            edgestrokec=col_edges,
            arrowlengthfrac=0.05,
        )
        return p
    else
        pos = collect(Iterators.product(1:height, 1:width))[:]
        locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

        p = GraphPlot.gplot(sample,
            locs_x,
            reverse(locs_y),
            nodelabel=LG.vertices(sample),
            nodefillc=col_nodes,
            edgestrokec=col_edges,
            arrowlengthfrac=0.05,
        )
        return p
    end
end

function color_cycles(
        graph::LG.SimpleDiGraph{T}
) where {T}

    edge_map = edgemap(graph)
    edge_idx = T[]
    nodes = Set{T}()
    for cycle in LG.simplecycles(graph)
        union!(nodes, cycle)
        for (x, y) in zip(cycle, circshift(cycle, -1))
            push!(edge_idx, edge_map[LG.Edge(x, y)])
        end
    end

    col_nodes = [colorant"turquoise" for _ in 1:LG.nv(graph)]
    col_nodes[collect(nodes)] .= colorant"red"

    col_edges = [colorant"lightgray" for _ in 1:LG.ne(graph)]
    col_edges[edge_idx] .= colorant"red"

    return col_nodes, col_edges
end

function plot(
    sfg::SinkFreeGraph,
    sample,
    width::Int=0,
    height::Int=0
)
    col_nodes = [LG.outdegree(sample, v) == 0 ? colorant"red" : colorant"turquoise"
                 for v in LG.vertices(sample)]
    if width == 0 || height == 0
        p = GraphPlot.gplot(sample,
            nodelabel=LG.vertices(sample),
            nodefillc=col_nodes,
            arrowlengthfrac=0.05,
        )
        return p
    else
        pos = collect(Iterators.product(1:height, 1:width))[:]
        locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

        p = GraphPlot.gplot(sample,
            locs_x,
            reverse(locs_y),
            nodelabel=LG.vertices(sample),
            nodefillc=col_nodes,
            arrowlengthfrac=0.05,
        )
        return p
    end
end

function color_sinks(
    graph::LG.SimpleDiGraph{T}
) where {T}

    nodes = T[]
    edge_idx = T[]
    edge_map = edgemap(graph)
    for v in sink_nodes(graph)
        push!(nodes, v)
        for w in LG.inneighbors(graph, v)
            push!(edge_idx, edge_map[LG.Edge(w, v)])
        end
    end

    col_nodes = [Colors.colorant"turquoise" for _ in 1:LG.nv(graph)]
    col_nodes[nodes] .= Colors.colorant"red"

    col_edges = [Colors.colorant"lightgray" for _ in 1:LG.ne(graph)]
    col_edges[edge_idx] .= Colors.colorant"red"

    return col_nodes, col_edges
end
