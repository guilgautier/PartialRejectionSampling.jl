using PartialRejectionSampling
const PRS = PartialRejectionSampling

β₀ = 0.1
r = 0.05  # interaction range = 2*radius
b = β₀ / (π * (r/2)^2)

c, w = [0.0, 0.0], 1.0
win = PRS.SquareWindow(c, w)

hc = PRS.HardCorePointProcess(b, r, win)

seed = -1
@time sample = PRS.generate_sample_prs(hc; rng=seed)
# @time sample = PRS.generate_sample_grid_prs(hc; rng=-1)

p = PRS.plot(hc, sample; title="")
