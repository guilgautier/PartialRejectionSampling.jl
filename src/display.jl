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
    #     nodelabel=LG.vertices(g),
    #     arrowlengthfrac=0.05
    #     edgestrokec=col_edges

    return p
end
