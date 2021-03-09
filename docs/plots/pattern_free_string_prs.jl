n = 100
alph = ["A", "C", "G", "T"]
pattern = "ACAC"
pfs = PRS.PatternFreeString(alph, pattern)

using Random
rng = Random.MersenneTwister(123)
@time PRS.generate_sample_prs(rng, pfs, n)
