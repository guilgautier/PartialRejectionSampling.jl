include(joinpath(@__DIR__, "pedagogy_graph.jl"))
using Random

dims = [5, 5]
g = LG.grid(dims)
β = 0.5

pp = PRS.HardCoreGraph(g, β)
rng = Random.MersenneTwister(123)

hard_core_graph_pedagogy(pp, dims, rng)
