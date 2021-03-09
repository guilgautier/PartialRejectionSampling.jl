# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_graph.jl")

dims = [5, 5]
g = LG.grid(dims)

sfg = PRS.SinkFreeGraph(g)

using Random
rng = Random.MersenneTwister(123)
@time sample = PRS.generate_sample_prs(rng, sfg)

p = plot(sfg, sample, dims, "docs/plots/output/sink_free_graph.pdf")
