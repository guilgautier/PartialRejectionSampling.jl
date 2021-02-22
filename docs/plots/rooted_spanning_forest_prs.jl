# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_graph.jl")

dims = [5, 5]
g = LG.grid(dims)
roots = [1, 5]

rsf = PRS.RootedSpanningForest(g, roots)

seed = 123
@time sample = PRS.generate_sample_prs(rsf; rng=seed)

p = plot(rsf, sample, dims...; file="docs/plots/output/rooted_spanning_forest.pdf")
