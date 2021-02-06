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
include("window.jl")

# sampling
include("dominated_cftp.jl")
include("grid_prs.jl")

# Spatial point processes
include("poisson.jl")
include("strauss.jl")
include("hard_core_spatial.jl")

# Graph point processes
include("ising.jl")
include("hard_core_graph.jl")
include("rooted_spanning_forest.jl")
include("sink_free_graph.jl")

# Misc
include("pattern_free_string.jl")

end
