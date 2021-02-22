# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_graph.jl")

dims = [5, 5]
g = LG.grid(dims)

sfg = PRS.SinkFreeGraph(g)

seed = 123
@time sample = PRS.generate_sample_prs(sfg; rng=seed)

p = plot(sfg, sample, dims...; file="docs/plots/output/sink_free_graph.pdf")
