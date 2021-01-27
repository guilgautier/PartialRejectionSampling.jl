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

p = PRS.plot(rsf, sample, width, height)
