include(joinpath(@__DIR__, "pedagogy_graph.jl"))
using Random

dims = [5, 5]
g = LG.grid(dims)
roots = [8]
pp = PRS.RootedSpanningForest(g, roots)

rng = Random.MersenneTwister(7)
rooted_spanning_tree_padagogy(pp, dims, rng)
