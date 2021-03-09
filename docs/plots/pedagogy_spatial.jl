using PartialRejectionSampling
using Distributions

using Plots

# Pedagogy
function pedagogy_generate_sample_prs(
    pp::PRS.HardCorePointProcess{T};
    win::Union{Nothing,PRS.AbstractWindow}=nothing,
    rng=-1
)::Vector{T} where {T}

    path(i) = joinpath("docs/plots/output/hard_core_spatial", join([lpad(i, 3, "0"), ".pdf"]))

    window_ = win === nothing ? PRS.window(pp) : win

    n = rand(rng, Distributions.Poisson(pp.β * PRS.volume(window_)))
    points = Matrix{Float64}(undef, PRS.dimension(pp), n)
    for x in eachcol(points)
        x .= rand(rng, window_)
    end

    i = 0
    while true

        p = pedagogy_plot(points, true, 0, "white", pp.window)
        Plots.savefig(p, path(i))

        p = pedagogy_plot(points, true, pp.r/2, "white", pp.window)
        i+=1
        Plots.savefig(p, path(i))

        bad = vec(any(PRS.pairwise_distances(points) .< pp.r, dims=2))
        !any(bad) && break

        p = pedagogy_plot(points[:, .!bad], true, 0, "white", pp.window)
        pedagogy_plot!(p, points[:, bad], false, pp.r/2, "red", pp.window)
        i+=1
        Plots.savefig(p, path(i))

        p = pedagogy_plot(points[:, .!bad], true, 0, "white", pp.window)
        pedagogy_plot!(p, points[:, bad], false, pp.r, "orange", pp.window)
        i+=1
        Plots.savefig(p, path(i))

        resampled = PRS.generate_sample_poisson_union_balls(rng, pp.β, points[:, bad], pp.r, window_)
        pedagogy_plot!(p, resampled, true, 0, "blue", pp.window)
        i+=1
        Plots.savefig(p, path(i))

        points = hcat(points[:, .!bad], resampled)
        i += 1
    end

    i += 1
    p = pedagogy_plot(points, true, 0, "white", pp.window)
    Plots.savefig(p, path(i))

    return [points[:, i] for i in 1:size(points, 2)]
end

function pedagogy_plot!(
    p,
    points,
    show_center=true,
    radius=0,
    color="white",
    window=SquareWindow(zeros(2), 1)
)
    for x in (points isa Matrix ? eachcol(points) : points)
        if radius > 0
            plot_disk!(p, x, radius, color)
        end
        if show_center
            Plots.scatter!(p, [x[1]], [x[2]], markersize=2, color=color, grid=false)
        end
    end

    Plots.xlims!(window.c[1], window.c[1] + window.w)
    Plots.ylims!(window.c[2], window.c[2] + window.w)

    return p
end

function pedagogy_plot(
    points,
    show_center=true,
    radius=0,
    color="white",
    window=SquareWindow(zeros(2), 1)
)

    p = Plots.plot([0], [0],
            label="", legend=false,
            color="white",
            linewidth=0,
            aspect_ratio=:equal,
            grid=false,
            title="")

    pedagogy_plot!(p, points, show_center, radius, color, window)

    return p
end

function plot_disk!(p, center, radius, color=:white)
    θ = LinRange(0, 2π, 15)
    x, y = center[1] .+ radius .* cos.(θ), center[2] .+ radius .* sin.(θ)

    Plots.plot!(p, x, y, seriestype=[:shape], c=color, linewidth=0.5, legend=false, fillapha=0.2)

    return p
end

function plot(
    pp::PRS.AbstractSpatialPointProcess,
    points;
    title=""
)
    p = Plots.plot([0], [0],
            label="", legend=false,
            color="white",
            linewidth=0,
            aspect_ratio=:equal,
            grid=false,
            title=title)

    θ = collect(range(0, 2π, length=15))
    rad = pp.r / 2  # radius = interaction range / 2
    circ_x, circ_y = rad .* cos.(θ), rad .* sin.(θ)

    for x in points
        Plots.plot!(x[1] .+ circ_x,
                    x[2] .+ circ_y,
                    color="black",
                    linewidth=0.5,
                    grid=false)
    end

    win = pp.window
    Plots.xlims!(win.c[1], win.c[1] + win.w[1])
    Plots.ylims!(win.c[2], win.c[2] + (win.w isa Number ? win.w[1] : win.w[2]))

    return p
end
