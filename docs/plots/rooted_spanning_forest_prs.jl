include(joinpath(@__DIR__, "pedagogy_graph.jl"))
using Random

dims = [5, 5]
g = LG.grid(dims)
roots = [1, 5]

rsf = PRS.RootedSpanningForest(g, roots)

rng = Random.MersenneTwister(123)
@time sample = PRS.generate_sample_prs(rng, rsf)

path = joinpath(@__DIR__, "output", "rooted_spanning_tree", "rooted_spanning_forest.pdf")
p = plot(rsf, sample, dims, path)
