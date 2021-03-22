using PartialRejectionSampling
using Random

n = 100
alphabet = ["A", "C", "G", "T"]
pattern = "ACAC"
pfs = PRS.PatternFreeString(pattern, alphabet)

rng = Random.MersenneTwister(123)
@time PRS.generate_sample_prs(rng, pfs, n)
