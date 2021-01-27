## Dominated Coupling From The Past (dCFTP)
#  - [Kendall \& Moller's orginal formulation of dCFTP](https://www.researchgate.net/publication/  2821877_Perfect_Metropolis-Hastings_simulation_of_locally_stable_point_processes)
#  - Huber, Perfect Simulation
#  - [Kendall's notes on perfect simulation](https://warwick.ac.uk/fac/sci/statistics/staff/  academic-research/kendall/personal/ppt/428.pdf)
#
# Requirement: the target spatial point process must have the following methods
#  - function papangelou_conditional_intensity end
#  - function upper_bound_papangelou_conditional_intensity end
#  - function window end

function generate_sample_dcftp(
        pp::AbstractSpatialPointProcess{T};
        n₀::Int=1,
        win::Union{Nothing,AbstractWindow}=nothing,
        rng=-1
)::Vector{T} where {T}

    @assert n₀ >= 1
    rng = getRNG(rng)

    window_ = win === nothing ? window(pp) : win
    β = upper_bound_papangelou_conditional_intensity(pp)
    birth_rate = β * volume(window_)

    # Dominating process
    hp = HomogeneousPoissonPointProcess(β, window_)
    D = Set{T}(eachcol(generate_sample(hp; rng=rng)))

    M = Float64[]  # Marking process
    R = T[]        # Recording process

    steps = -1:-1:-n₀
    while true
        backward_update!(D, M, R, steps, birth_rate, window_; rng=rng)
        coupling, L = forward_coupling(D, M, R, pp, β)
        coupling && return collect(L)
        steps = (steps.stop-1):-1:(2*steps.stop)
    end
end

function backward_update!(
        D::Set{T},          # Dominating process
        M::Vector{Float64}, # Marking process
        R::Vector{T},       # Recording process
        steps::StepRange,   # Number of backward steps
        birth_rate::Real,
        win::AbstractWindow;
        rng=-1
) where {T}
    rng = getRNG(rng)
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
            x = rand(win; rng=rng)
            push!(D, x)
            pushfirst!(R, x)
            pushfirst!(M, 0.0)
        end
    end
end

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
