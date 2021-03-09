# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_graph.jl")

dims = [5, 5]
g = LG.grid(dims)
roots = [1, 5]

rsf = PRS.RootedSpanningForest(g, roots)

using Random
rng = Random.MersenneTwister(123)
@time sample = PRS.generate_sample_prs(rng, rsf)

p = plot(rsf, sample, dims, "docs/plots/output/rooted_spanning_forest.pdf")
