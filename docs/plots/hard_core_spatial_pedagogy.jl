# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_spatial.jl")

β₀ = 0.3
r = 0.1  # interaction range = 2*radius
b = β₀ / (π * (r/2)^2)

c, w = [0.0, 0.0], 1.0
win = PRS.SquareWindow(c, w)

hc = PRS.HardCorePointProcess(b, r, win)

seed = 123
pedagogy_generate_sample_prs(hc; rng=seed)
