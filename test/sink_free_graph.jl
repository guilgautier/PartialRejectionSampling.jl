using PartialRejectionSampling
const PRS = PartialRejectionSampling

using LightGraphs
const LG = LightGraphs

width, height = 5, 5
g = LG.grid([width, height])

sfg = PRS.SinkFreeGraph(g)
seed = -1
@time sample = PRS.generate_sample_prs(sfg; rng=seed)

using GraphPlot
using Colors

pos = collect(Iterators.product(1:height, 1:width))[:]
locs_x, locs_y = map(x->x[1], pos), map(x->x[2], pos)

sample = PRS.generate_sample_prs(sfg; rng=seed)
cols = [LG.outdegree(sample, v) == 0 ? colorant"red" : colorant"turquoise"
        for v in LG.vertices(sample)]
p = GraphPlot.gplot(sample,
    locs_x, locs_y,
    nodelabel=LG.vertices(sample),
    nodefillc=cols,
    arrowlengthfrac=0.05
)
display(p)
