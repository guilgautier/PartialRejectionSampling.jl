struct Ising{T<:Integer} <: AbstractGraphPointProcess{T}
    "Grid dimension"
    dims::Vector{T}
    "Grid graph"
    g::LG.SimpleGraph{T}
    "Magnetization"
    h::Union{Float64, Vector{Float64}}
    "Interaction"
    J::Float64
end

function Ising(
        dims::Vector{T},
        periodic::Bool,
        h::Union{Real, Vector{Real}},
        J::Real
) where {T<:Integer}
    h isa AbstractVector && @assert length(h) == prod(dims)
    g = LG.grid(dims; periodic=periodic)
    return Ising(dims, g, float.(h), float(J))
end

function Ising(
        dims::Vector{T},
        periodic::Bool,
        J::Real
) where {T<:Integer}
    return Ising(dims, periodic, 0, J)
end

## Sampling

include("ising/ising_sampling_gibbs_perfect.jl")
include("ising/ising_sampling_grid_prs.jl")
