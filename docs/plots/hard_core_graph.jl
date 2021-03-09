# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_graph.jl")

dims = [5, 5]
g = LG.grid(dims)
β = 1.0

hcg = PRS.HardCoreGraph(g, β)

using Random
rng = Random.MersenneTwister(123)
@time sample = PRS.generate_sample_prs(rng, hcg)

p = plot(hcg, sample, dims, "docs/plots/output/hard_core_graph.pdf")

@time sample = PRS.generate_sample_prs(hcg)
