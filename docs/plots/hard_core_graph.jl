include(joinpath(@__DIR__, "pedagogy_graph.jl"))
using Random

dims = [5, 5]
g = LG.grid(dims)
β = 1.0

hcg = PRS.HardCoreGraph(g, β)

rng = Random.MersenneTwister(123)
@time sample = PRS.generate_sample_prs(rng, hcg)

path = joinpath(@__DIR__, "output", "hard_core_graph", "hard_core_graph.pdf")
p = plot(hcg, sample, dims, path)

# @time sample = PRS.generate_sample_prs(hcg)
