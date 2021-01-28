using PartialRejectionSampling
const PRS = PartialRejectionSampling

β₀ = 0.1
r = 0.05  # interaction range = 2*radius
β, γ = β₀ / (π * (r/2)^2), 0.1  # γ = 0 ≡ Hard core

c, w = [0.0, 0.0], 2.0
win = PRS.SquareWindow(c, w)

strauss = PRS.StraussPointProcess(β, γ, r, win)

seed = -1
@time sample = PRS.generate_sample_grid_prs(strauss; rng=seed)

p = PRS.plot(strauss, sample)
