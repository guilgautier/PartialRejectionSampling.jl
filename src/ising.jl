struct Ising{T<:Integer} <: AbstractGraphPointProcess{T}
    "Grid dimension"
    dims::Vector{T}
    "Grid graph"
    g::LG.SimpleGraph{T}
    "Magnetization"
    h::Union{Float64,Vector{Float64}}
    "Interaction"
    J::Float64
end

function Ising(
        dims::Vector{T},
        periodic::Bool,
        J::Real,
        h::Real=0,
) where {T<:Integer}
    g = LG.grid(dims; periodic=periodic)
    return Ising{T}(dims, g, float(h), J)
end

function Ising(
        dims::Vector{T},
        periodic::Bool,
        J::Real,
        h::Vector{Real},
) where {T<:Integer}
    @assert length(h) == prod(dims)
    g = LG.grid(dims; periodic=periodic)
    return Ising{T}(dims, g, float.(h), J)
end

## Sampling

include("ising/ising_sampling_gibbs_perfect.jl")
include("ising/ising_sampling_grid_prs.jl")
