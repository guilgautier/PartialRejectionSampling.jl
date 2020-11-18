# using Revise

import PartialRejectionSampling
const PRS = PartialRejectionSampling

dims = [5, 5] # if > (14, 14) the display becomes all black, don't know why !
periodic = false
H, J = 0.0, 0.001

ising = PRS.Ising(dims, periodic, H, J)

println("\nPartial Rejection Sampling\n")

seed = -1

state, cnt = PRS.prs(ising; rng=seed)
energy = PRS.energy(ising, state)
println("Number of resampling steps:\n$cnt")
println("Energy:\n$energy")
PRS.plot(ising, state)

println("\nPerfect Gibbs sampler\n")

seed = -1

state, cnt = PRS.gibbs_perfect(ising; rng=seed)
energy = PRS.energy(ising, state)
println("Number of resampling steps:\n$cnt")
println("Energy:\n$energy")
PRS.plot(ising, state)

