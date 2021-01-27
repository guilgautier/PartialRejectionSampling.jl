abstract type AbstractPointProcess{T} end
Base.eltype(pp::AbstractPointProcess{T}) where {T} = T

abstract type AbstractSpatialPointProcess{T<:Vector{Float64}} <: AbstractPointProcess{T} end

window(pp::AbstractSpatialPointProcess) = pp.window
dimension(pp::AbstractSpatialPointProcess) = dimension(window(pp))

abstract type AbstractGraphPointProcess{T} <: AbstractPointProcess{T} end

function generate_sample end
