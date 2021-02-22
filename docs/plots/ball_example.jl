using PartialRejectionSampling
const PRS = PartialRejectionSampling

c, w = [0.0, 0.0], 10.0
win = PRS.SquareWindow(c, w)

r = 4
hc = PRS.HardCorePointProcess(b, r, win)

@time pts = rand(PRS.BallWindow(c, r), 10)

using Plots

Plots.scatter(pts[1, :], pts[2, :], aspect_ratio=:equal)
