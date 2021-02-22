# const PRS = PartialRejectionSampling, Plots, Colors, LS = LazySets
includet("docs/plots/pedagogy_spatial.jl")

β₀ = 0.1
r = 0.05  # interaction range = 2*radius
β, γ = β₀ / (π * (r/2)^2), 0.1  # γ = 0 ≡ Hard core

β₀ = 0.3
r = 0.07  # interaction range = 2*radius
β, γ = β₀ / (π * (r/2)^2), 0.1  # γ = 0 ≡ Hard core

c, w = [0.0, 0.0], 1.0
win = PRS.SquareWindow(c, w)

strauss = PRS.StraussPointProcess(β, γ, r, win)

seed = 123
@time sample = PRS.generate_sample_grid_prs(strauss; rng=seed)

p = pedagogy_plot(sample,
    true,
    strauss.r/2,
    "white",
    win
)

Plots.savefig(p, "docs/plots/output/strauss/strauss_grid_prs.pdf")
