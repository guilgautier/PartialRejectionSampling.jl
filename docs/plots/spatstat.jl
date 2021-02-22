using PartialRejectionSampling;
const PRS = PartialRejectionSampling;

using RCall;
@rimport spatstat as rspatstat;

# Create the point process
β₀ = 0.1
r = 0.02  # interaction range = 2*radius
b = β₀ / (π * (r/2)^2)

c, w = zeros(2), 1.0
win = PRS.SquareWindow(c, w)

hc = PRS.HardCorePointProcess(b, r, win)

# Time comparison PRS / spatstat
@time PRS.generate_sample_prs(hc; rng=-1);
@time rspatstat.rHardcore(beta=hc.β, R=hc.r);

# K-plot with spatstat
R"""
library(spatstat)
sampl_spatstat = rHardcore(beta=$b, R=$r)
K <- Kest(sampl_spatstat)
plot(K)
"""

# K-plot with PRS
sampl_PRS = hcat(PRS.generate_sample_prs(hc; rng=-1)...)

R"""
library(spatstat)

tmp = $(hc.window.c) + $(hc.window.w)
xrange <- c($(hc.window.c[1]), tmp[1])
yrange <- c($(hc.window.c[2]), tmp[2])
window <- owin(xrange, yrange)

sampl <- ppp($(sampl_PRS[1, :]), $(sampl_PRS[2, :]), window)
K <- Kest(sampl)
plot(K)
"""
