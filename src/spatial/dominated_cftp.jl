"""
Implementation of dominated Coupling From The Past (dCFTP) developed by [KeMo99](@cite) and [KeMo00](@cite) for [`PRS.AbstractSpatialPointProcess`](@ref).

**See also**

- [Hub16](@cite)
- [Kendall's notes on perfect simulation](https://warwick.ac.uk/fac/sci/statistics/staff/  academic-research/kendall/personal/ppt/428.pdf)
"""

@doc raw"""
    papangelou_conditional_intensity(pp::AbstractSpatialPointProcess, x, X)

Compute the [Papangelou conditional intensity](https://en.wikipedia.org/wiki/Point_process#Papangelou_intensity_function) of `pp`, as the ratio of densities ``\frac{f(X \cup x)}{f(X)}``.

**See also**

- Section 6.1.1 [MoWa04](@cite)
"""
function papangelou_conditional_intensity(pp::AbstractSpatialPointProcess, x, X) end

"""
    upper_bound_papangelou_conditional_intensity(pp::AbstractSpatialPointProcess)

Compute an upper bound on the [`PRS.papangelou_conditional_intensity`](@ref) of `pp`

**See also**

- Equation 2.1 [KeMo99](@cite)
"""
function upper_bound_papangelou_conditional_intensity(pp::AbstractSpatialPointProcess) end

"""
    isrepulsive(pp::AbstractSpatialPointProcess)

Property of a [`PRS.AbstractSpatialPointProcess`](@ref).
"""
function isrepulsive(pp::AbstractSpatialPointProcess) end

"""
    isattractive(pp::AbstractSpatialPointProcess)

Property of a [`PRS.AbstractSpatialPointProcess`](@ref).
"""
function isattractive(pp::AbstractSpatialPointProcess) end

"""
    generate_sample_dcftp(
        [rng::Random.AbstractRNG,]
        pp::AbstractSpatialPointProcess{T};
        win::Union{Nothing,AbstractSpatialWindow}=nothing,
        n₀::Int=1
    )::Vector{T} where {T}

Generate an exact sample from a spatial point process `pp` on window `win` using dominated coupling from the past.

- Default window (`win=nothing`) is `window(pp)=pp.window`
- Initial coalescence check performed after `n₀` steps.
- Seed or random number generator is passed via `rng`.

**See also**

- [KeMo99](@cite), [KeMo00](@cite)
- Section 11.2.6 [MoWa04](@cite)
"""
function generate_sample_dcftp(
    rng::Random.AbstractRNG,
    pp::AbstractSpatialPointProcess{T},
    win::Union{Nothing,AbstractSpatialWindow}=nothing,
    n₀::Int=1
)::Vector{T} where {T}
    @assert n₀ >= 1

    window_ = isnothing(win) ? window(pp) : win
    β = upper_bound_papangelou_conditional_intensity(pp)
    birth_rate = β * volume(window_)

    # Dominating process
    hppp = HomogeneousPoissonPointProcess(β, window_)
    D = Set{T}(eachcol(generate_sample(rng, hppp)))

    M = Float64[]  # Marking process
    R = T[]        # Recording process

    steps = -1:-1:-n₀
    while true
        backward_extend!(rng, D, M, R, steps, birth_rate, window_)
        coupling, L = forward_coupling(D, M, R, pp, β)
        coupling && return collect(L)
        steps = (steps.stop-1):-1:(2*steps.stop)
    end
end

function generate_sample_dcftp(pp::AbstractPointProcess, args...)
    return generate_sample_dcftp(Random.default_rng(), pp, args...)
end

"""
    backward_extend!(
        rng::Random.AbstractRNG,
        D::Set{T},          # Dominating process
        M::Vector{Float64}, # Marking process
        R::Vector{T},       # Recording process
        steps::StepRange,   # Number of backward steps
        birth_rate::Real,
        win::AbstractSpatialWindow
    ) where {T}

Sample from the dominating birth-death process backwards in time according to `steps`.
The marks `M` and the points `R` which were added (uniform mark) / deleted (mark=0) along the run are recorded (`pushfirst!`).
The final state of the dominating process is updated to `D`.

**See also**

- `SBDevolve()` [KeMo99](@cite)
"""
function backward_extend!(
    rng::Random.AbstractRNG,
    D::Set{T},          # Dominating process
    M::Vector{Float64}, # Marking process
    R::Vector{T},       # Recording process
    steps::StepRange,   # Number of backward steps
    birth_rate::Real,
    win::AbstractSpatialWindow
) where {T}
    for _ in steps
        card_D = length(D)
        if rand(rng) < card_D / (birth_rate + card_D)
            # forward death (delete) ≡ backward birth (pushfirst)
            x = rand(rng, D)
            delete!(D, x)
            pushfirst!(R, x)
            pushfirst!(M, rand(rng))
        else
            # forward birth (push) ≡ backward death (pushfirst)
            x = rand(rng, win)
            push!(D, x)
            pushfirst!(R, x)
            pushfirst!(M, 0.0)
        end
    end
end

"""
    forward_coupling(
        D::Set{T},           # Dominating process
        M::Vector{Float64},  # Marking process
        R::Vector{T},        # Recording process
        pp::AbstractSpatialPointProcess{T},
        β::Real              # Upper bound on papangelou conditional intensity
    ) where {T}

Check if coalescence occured between the lower and upper bounding processes and return the state of the lower bounding process at time 0.

**See also**

- `SBDadd()` [KeMo99](@cite)
"""
function forward_coupling(
    D::Set{T},           # Dominating process
    M::Vector{Float64},  # Marking process
    R::Vector{T},        # Recording process
    pp::AbstractSpatialPointProcess{T},
    β::Real              # Upper bound on papangelou conditional intensity
) where {T}
    # L ⊆ X ⊆ U ⊆ D, where X is the target process
    L, U = empty(D), copy(D)
    for (m, x) in zip(M, R)
        if m > 0  # if birth occured in D
            if isrepulsive(pp)
                if m < papangelou_conditional_intensity(pp, x, U) / β
                    push!(L, x)
                    push!(U, x)
                elseif m < papangelou_conditional_intensity(pp, x, L) / β
                    push!(U, x)
                end
            elseif isattractive(pp)
                if m < papangelou_conditional_intensity(pp, x, L) / β
                    push!(L, x)
                    push!(U, x)
                elseif m < papangelou_conditional_intensity(pp, x, U) / β
                    push!(U, x)
                end
            else
                error("The current implementation of Dominated Coupling From The Past "
                       * "requires the target point process to be attractive or repulsive.")
            end
        else  # if death occured in D
            delete!(L, x)
            delete!(U, x)
        end
    end
    # Check coalescence L = U (Since L ⊆ U: L = U ⟺ |L| = |U|)
    return length(L) == length(U), L
end
