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
        state
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
    return p
end

function plot(
        hcg::HardCoreGraph,
        state,
        width::Int=0,
        height::Int=0
)
    col_nodes = [Colors.colorant"turquoise" for _ in 1:LG.nv(hcg.g)]
    col_nodes[state] .= Colors.colorant"red"

    if width == 0 || height == 0
        p = GraphPlot.gplot(hcg.g,
                nodelabel=LG.vertices(hcg.g),
                nodefillc=col_nodes
                )
        return p
    else
        pos = collect(Iterators.product(1:height, 1:width))[:]
        locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

        p = GraphPlot.gplot(hcg.g,
                locs_x, locs_y,
                nodelabel=LG.vertices(hcg.g),
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
            locs_x, locs_y,
            nodelabel=LG.vertices(sample),
            nodefillc=col_nodes,
            edgestrokec=col_edges,
            arrowlengthfrac=0.05,
        )
        return p
    end
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
            locs_x, locs_y,
            nodelabel=LG.vertices(sample),
            nodefillc=col_nodes,
            arrowlengthfrac=0.05,
        )
        return p
    end
end
