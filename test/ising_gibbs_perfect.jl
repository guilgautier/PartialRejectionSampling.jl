import PartialRejectionSampling
const PRS = PartialRejectionSampling

dims = [15, 15] # if > (14, 14) the display becomes all black, don't know why !
periodic = false
H, J = 0.0, 0.01

seed = -1

ising = PRS.Ising(dims, periodic, H, J)

@time sample = PRS.generate_sample_gibbs_perfect(ising; rng=seed)

PRS.plot(ising, sample)
