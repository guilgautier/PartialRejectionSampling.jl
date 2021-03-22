include(joinpath(@__DIR__, "pedagogy_graph.jl"))
using Random

dims = [5, 5]
g = LG.grid(dims)

pp = PRS.SinkFreeGraph(g)
rng = Random.MersenneTwister(7)
sink_free_graph_pedagogy(pp, dims, rng)
