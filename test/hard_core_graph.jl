using PartialRejectionSampling
const PRS = PartialRejectionSampling

using LightGraphs
const LG = LightGraphs

width, height = 5, 5
g = LG.grid([width, height])
β = 1.0

hcg = PRS.HardCoreGraph(g, β)
seed = -1
@time sample = PRS.generate_sample_prs(hcg; rng=seed)

using GraphPlot
using Colors

pos = collect(Iterators.product(1:height, 1:width))[:]
locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

col_nodes = [Colors.colorant"turquoise" for _ in 1:LG.nv(g)]
col_nodes[sample] .= Colors.colorant"red"

p = gplot(g,
    locs_x, locs_y,
    nodelabel=LG.vertices(g),
    nodefillc=col_nodes,
#     arrowlengthfrac=0.05
#     edgestrokec=col_edges
)
display(p)
