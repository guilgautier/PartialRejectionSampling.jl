n = 100
alph = ["A", "C", "G", "T"]
pattern = "ACAC"
pfs = PRS.PatternFreeString(alph, pattern)
seed = 123
@time PRS.generate_sample_prs(pfs, n; rng=seed)
