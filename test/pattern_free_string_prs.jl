using PartialRejectionSampling
const PRS = PartialRejectionSampling

n = 1000
alph = ["A", "C", "G", "T"]
pattern = "ACAC"

seed = -1
@time PRS.generate_pattern_free_string(n, alph, pattern; rng=seed)
