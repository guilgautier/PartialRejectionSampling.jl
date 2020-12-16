using PartialRejectionSampling
const PRS = PartialRejectionSampling

using LightGraphs
const LG = LightGraphs

width, height = 5, 5
g = LG.grid([width, height])
roots = [1, 5]

rsf = PRS.RootedSpanningForest(g, roots)
seed = -1
@time sample = PRS.generate_sample_prs(rsf; rng=seed)

# Display

using GraphPlot
using Colors

pos = collect(Iterators.product(1:height, 1:width))[:]
locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

col_nodes, col_edges = PRS.color_cycles(sample)
col_nodes[roots] .= Colors.colorant"orange"

p = GraphPlot.gplot(sample,
    locs_x, locs_y,
    nodelabel=LG.vertices(sample),
    nodefillc=col_nodes,
    edgestrokec=col_edges,
    arrowlengthfrac=0.05,
)
display(p)
