include(joinpath(@__DIR__, "pedagogy_graph.jl"))
using Random

dims = [5, 5] # if > (14, 14) the display becomes all black, don't know why !
periodic = false
H, J = 0.0, 0.01

ising = PRS.Ising(dims, J, H; periodic=periodic)

using Random
rng = Random.MersenneTwister(123)
@time sample = PRS.generate_sample_gibbs_perfect(rng, ising)

path = joinpath(@__DIR__, "output", "ising", "ising_grid_prs.pdf")
p = plot(ising, sample, dims, path; nodelabel=nothing)
