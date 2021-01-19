module PartialRejectionSampling

using LinearAlgebra
const LA = LinearAlgebra

using Random, Distributions

using LightGraphs
const LG = LightGraphs

using SimpleWeightedGraphs
const SWG = SimpleWeightedGraphs

using Plots, GraphPlot, Colors

using Distances

using LazySets
const LS = LazySets

abstract type AbstractPointProcess{T} end
Base.eltype(pp::AbstractPointProcess{T}) where {T} = T
abstract type AbstractSpatialPointProcess{T<:Vector} <: AbstractPointProcess{T} end
abstract type AbstractGraphPointProcess{T} <: AbstractPointProcess{T} end

function generate_sample end

include("utils.jl")
include("window.jl")

# sampling
include("dominated_cftp.jl")
include("grid_prs.jl")

# Spatial point processes
include("strauss.jl")
include("hard_core_spatial.jl")

include("hard_core_graph.jl")
include("ising.jl")

# Graph point processes
include("rooted_spanning_forest.jl")
include("sink_free_graph.jl")

# Misc
include("pattern_free_string.jl")

# Display
include("display.jl")

end
