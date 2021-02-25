# const PRS = PartialRejectionSampling, Plots, GraphPlot, Colors, LS = LazySets, LG = LightGraphs
includet("docs/plots/pedagogy_graph.jl")

dims = [5, 5] # if > (14, 14) the display becomes all black, don't know why !
periodic = false
H, J = 0.0, 0.01

ising = PRS.Ising(dims, J, H; periodic=periodic)

seed = 123
@time sample = PRS.generate_sample_grid_prs(ising; rng=seed)

p = plot(ising, sample, dims...; file="docs/plots/output/ising/ising_grid_prs.pdf")