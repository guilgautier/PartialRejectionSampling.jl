include(joinpath(@__DIR__, "pedagogy_graph.jl"))
using Random

dims = [5, 5]
g = LG.grid(dims)

sfg = PRS.SinkFreeGraph(g)

rng = Random.MersenneTwister(123)
@time sample = PRS.generate_sample_prs(rng, sfg)

path = joinpath(@__DIR__, "output" "sink_free_graph", "sink_free_graph.pdf")
p = plot(sfg, sample, dims, path)
