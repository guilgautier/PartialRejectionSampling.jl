"""
    PatternFreeString{T<:String} <: AbstractPointProcess{T}

Container with fields `alphabet` and `pattern` used to generate a string made of characters from `alphabet` avoiding the prescribed `pattern`.
"""
struct PatternFreeString{T<:String} <: AbstractPointProcess{T}
    alphabet::Vector{T}
    pattern::T
end

function Base.show(io::IO, pp::PatternFreeString{T}) where {T}
    print(io, "PatternFreeString{$T}\n- alphabet = $(pp.alphabet)\n- pattern = $(pp.pattern)")
end

"""
    PatternFreeString(alphabet::Vector{String}, pattern::String)

Construct a [`PRS.PatternFreeString`](@ref).

```jldoctest; output = true
using PartialRejectionSampling
PRS.PatternFreeString(["A", "C", "G", "T"], "ATGTA")

# output

PatternFreeString{String}
- alphabet = ["A", "C", "G", "T"]
- pattern = ATGTA
```
"""
function PatternFreeString(alphabet::Vector{String}, pattern::String)
    @assert !isempty(alphabet)
    @assert !isempty(pattern)
    for p in split(pattern, "")
        p ∉ alphabet && throw(DomainError(p, "pattern is not made of characters from alphabet"))
    end
    return PatternFreeString{String}(alphabet, pattern)
end

"""
    generate_sample(
        [rng=Random.AbstractRNG,]
        pp::PatternFreeString{T},
        size::Int
    )::T where {T<:String}

Generate a string uniformly at random among all strings made of characters from `pp.alphabet` with no occurence of the pattern `pp.pattern`.

Default sampler is [`PRS.generate_sample_prs`](@ref).
"""
function generate_sample(
    rng::Random.AbstractRNG,
    pp::PatternFreeString{T},
    size::Int
)::T where {T}
    return generate_sample_prs(rng, pp, size)
end

# function generate_sample(pp::PatternFreeString, size::Int)
#     return generate_sample_prs(Random.default_rng(), pp, size)
# end

"""
    generate_sample_prs(
        [rng::Random.AbstractRNG,]
        pp::PatternFreeString{T},
        size::Int
    )::T where {T<:String}

Generate a string uniformly at random among all strings made of characters from `pp.alphabet` with no occurence of the pattern `pp.pattern`, using a tailored version of Partial Rejection Sampling (PRS) derived by [GiAmWe18](@cite).

```@example
using PartialRejectionSampling
pp = PRS.PatternFreeString(["A", "C", "G", "T"], "ATGTA")
PRS.generate_sample_prs(pp, 20)
```
"""
function generate_sample_prs(
    rng::Random.AbstractRNG,
    pp::PatternFreeString{T},
    size::Int
)::T where {T<:String}
    @assert size > 0
    return _generate_sample_pattern_free_string_prs(rng, pp.alphabet, pp.pattern, size)
end

# function generate_sample_prs(
#     pp::PatternFreeString,
#     size::Int
# )
#     return generate_sample(Random.default_rng(), pp, size)
# end

"""
    _generate_pattern_free_string_prs(
        rng::Random.AbstractRNG,
        alphabet::Vector{T},
        pattern::T,
        size::Int
    )::T where {T<:AbstractString}

Generate a string uniformly at random among all strings made of characters from `alphabet` with no occurence of the pattern `pattern`, using a tailored version of Partial Rejection Sampling (PRS) derived by [GiAmWe18](@cite)
"""
function _generate_sample_pattern_free_string_prs(
    rng::Random.AbstractRNG,
    alphabet::Vector{T},
    pattern::T,
    size::Int
)::T where {T<:AbstractString}

    @assert size > 0
    @assert !isempty(alphabet)
    @assert !isempty(pattern)

    pref_suff = find_prefix_suffix(pattern)

    pp = fill("", size)
    resample_indices = Set(1:size)

    while !isempty(resample_indices)
        generate_sample!(rng, pp, resample_indices, alphabet)
        resample_indices = find_characters_to_resample(pp, pattern, pref_suff)
    end
    return join(pp)
end

find_prefix_suffix(s::String) = [i for i in 1:div(length(s), 2) if s[1:i] == s[end-i+1:end]]

"""
    generate_sample!(
        [rng::Random.AbstractRNG,]
        string_vec::Vector{T},
        indices,
        alphabet::Vector{T}
    ) where {T<:AbstractString}

Generate a character uniformly at random from `alphabet` at positions prescribed by `indices` in `string_vec`.
"""
function generate_sample!(
    rng::Random.AbstractRNG,
    string_vec::Vector{T},
    indices,
    alphabet::Vector{T}
) where {T<:AbstractString}
    for i in indices
        string_vec[i] = rand(rng, alphabet)
    end
end

"""
    find_bad_ranges(
        pattern::T,
        string::T
    )::Vector{UnitRange} where {T<:AbstractString}

Identify where `pattern` occur in `string` and return the corresponding ranges of indices.
"""
function find_bad_ranges(
    pattern::T,
    string::T
)::Vector{UnitRange} where {T<:AbstractString}

    bad_ranges = UnitRange{Int64}[]

    matches = eachmatch(Regex(pattern), string, overlap=true)
    isempty(matches) && return bad_ranges

    p = length(pattern)
    m, _ = iterate(matches)
    x1, y1 = m.offset, m.offset + p - 1

    for m in Iterators.drop(matches, 1)
        x2, y2 = m.offset, m.offset + p - 1
        if x2 <= y1 + 1
            y1 = y2
        else
            push!(bad_ranges, x1:y1)
            x1, y1 = x2, y2
        end
    end
    push!(bad_ranges, x1:y1)
    return bad_ranges
end

"""
    find_characters_to_resample(
        string_vec::Vector{T},
        pattern::T,
        pref_suff::Vector{U}
    )::Vector{U} where {T<:String, U<:Int}

Identify the set of events to be resampled as constructed by Algorithm 5 in [GuJeLi19](@cite) as part of the Partial Rejection Sampling (PRS) method.
Return the indices of the variables involved in the corresponding events.

**See also**

- [`PRS.find_bad_ranges`](@ref)
"""
function find_characters_to_resample(
    string_vec::Vector{T},
    pattern::T,
    pref_suff::Vector{U}
)::Vector{U} where {T<:String, U<:Int}

    # Extremal case
    isempty(pref_suff) && return vcat(findall(pattern, join(string_vec), overlap=false)...)

    # General case
    bad_ranges = find_bad_ranges(pattern, join(string_vec))
    isempty(bad_ranges) && return vcat(bad_ranges...)

    p, s = length(pattern), length(string_vec)
    tmp = fill("", s)

    for bad_range in bad_ranges
        tmp[bad_range] = string_vec[bad_range]
        start, stop = bad_range.start, bad_range.stop
        for ps in pref_suff
            flag_left = flag_right = false
            if !flag_left
                I = (start - p + ps):(start - 1)
                if I.start >= 1
                    flag_left = startswith(pattern, join(tmp[I]))
                    if flag_left
                        tmp[I] = string_vec[I]
                    end
                end
            end
            if !flag_right
                J = (stop + 1):(stop + p - ps)
                if J.stop <= s
                    flag_right = endswith(pattern, join(tmp[J]))
                    if flag_right
                        tmp[J] = string_vec[J]
                    end
                end
            end
            flag_left && flag_right && break
        end
    end
    return findall(!isempty, tmp)
end
