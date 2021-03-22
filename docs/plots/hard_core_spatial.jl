include(joinpath(@__DIR__, "pedagogy_spatial.jl"))
using Random

β₀ = 0.1
r = 0.05  # interaction range = 2*radius
b = β₀ / (π * (r/2)^2)

c, w = [0.0, 0.0], 1.0
win = PRS.SquareWindow(c, w)

hc = PRS.HardCorePointProcess(b, r, win)

rng = Random.MersenneTwister(123)
@time sample = PRS.generate_sample_prs(rng, hc)
# @time sample = PRS.generate_sample_grid_prs(rng, hc)

p = plot(hc, sample; title="")
