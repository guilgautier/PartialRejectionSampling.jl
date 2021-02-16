module PartialRejectionSampling

# Imports

using LinearAlgebra
const LA = LinearAlgebra

using Random, Distributions

using LightGraphs
const LG = LightGraphs

using SimpleWeightedGraphs
const SWG = SimpleWeightedGraphs

using Distances

using SpecialFunctions
const SF = SpecialFunctions

# Exports

const PRS = PartialRejectionSampling
export PRS

export AbstractWindow
export AbstractSpatialWindow,
    RectangleWindow,
    SquareWindow,
    BallWindow
export AbstractDiscreteWindow,
    GraphNode

export AbstractPointProcess
export AbstractSpatialPointProcess,
    HomogeneousPoissonPointProcess,
    HardCorePointProcess,
    StraussPointProcess
export AbstractGraphPointProcess,
    Ising,
    HardCoreGraph,
    RootedSpanningForest,
    SinkFreeGraph
export PatternFreeString

export generate_sample,
    generate_sample_prs,
    generate_sample_grid_prs,
    generate_sample_dcftp,         # For Spatial point processes
    generate_sample_gibbs_perfect  # For the Ising model

# Code inclusions

include("common.jl")
include("utils.jl")

include(joinpath("spatial", "window.jl"))
include(joinpath("graph", "window.jl"))

include("grid_prs.jl")

## Spatial point processes
include(joinpath("spatial", "dominated_cftp.jl"))
include(joinpath("spatial", "poisson.jl"))
include(joinpath("spatial", "hard_core.jl"))
include(joinpath("spatial", "strauss.jl"))

## Graph point processes
include(joinpath("graph", "ising.jl"))
include(joinpath("graph", "hard_core.jl"))
include(joinpath("graph", "rooted_spanning_forest.jl"))
include(joinpath("graph", "sink_free_graph.jl"))

# Misc
include(joinpath("misc", "pattern_free_string.jl"))

end
