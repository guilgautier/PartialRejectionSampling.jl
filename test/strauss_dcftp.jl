using PartialRejectionSampling
const PRS = PartialRejectionSampling

# Experimental setup refers to slide 9 of [E. Rubak](https://www-ljk.imag.fr/membres/Jean-Francois.Coeurjolly/documents/lecture4.pdf)
β, γ, r = 2.0, 0.2, 0.7  # γ = 0 ≡ Hard core

c, w = [0.0, 0.0], 10.0
win = PRS.SquareWindow(c, w)

strauss = PRS.StraussPointProcess(β, γ, r, win)

seed = -1
@time sample = PRS.generate_sample_dcftp(strauss; rng=seed)

p = PRS.plot(strauss, sample; title="")
