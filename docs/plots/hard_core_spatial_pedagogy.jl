include(joinpath(@__DIR__, "pedagogy_spatial.jl"))
using Random

β₀ = 0.3
r = 0.1  # interaction range = 2*radius
b = β₀ / (π * (r/2)^2)

c, w = [0.0, 0.0], 1.0
win = PRS.SquareWindow(c, w)

hc = PRS.HardCorePointProcess(b, r, win)

rng = Random.MersenneTwister(123)
hard_core_spatial_pedagogy(hc; rng=rng)
