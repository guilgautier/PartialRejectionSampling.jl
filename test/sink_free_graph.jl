using PartialRejectionSampling
const PRS = PartialRejectionSampling

using LightGraphs
const LG = LightGraphs

width, height = 5, 5
g = LG.grid([width, height])

sfg = PRS.SinkFreeGraph(g)

seed = -1
@time sample = PRS.generate_sample_prs(sfg; rng=seed)

p = PRS.plot(sfg, sample, width, height)
