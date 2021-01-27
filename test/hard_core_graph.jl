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

p = PRS.plot(hcg, sample, width, height)
