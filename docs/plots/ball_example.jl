using PartialRejectionSampling
const PRS = PartialRejectionSampling

c, r = rand(2), 10

using Random
rng = Random.MersenneTwister(123)
@time pts = rand(rng, PRS.BallWindow(c, r), 10)
@time pts = rand(PRS.BallWindow(c, r), 10)

using Plots

Plots.scatter(pts[1, :], pts[2, :], aspect_ratio=:equal)
