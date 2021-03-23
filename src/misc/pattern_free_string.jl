"""
    PatternFreeString{T<:String} <: AbstractPointProcess{T}

Struct with fields
- `alphabet::Vector{String}`,
- `pattern::String`,

used to generate strings made of characters from `alphabet` avoiding the prescribed `pattern`.
"""
struct PatternFreeString{T<:String} <: AbstractPointProcess{T}
    pattern::T
    alphabet::Vector{T}
end

function Base.show(io::IO, pp::PatternFreeString{T}) where {T}
    print(io, "PatternFreeString{$T}\n- pattern = $(pp.pattern)\n- alphabet = $(pp.alphabet)")
end

"""
    PatternFreeString(alphabet::Vector{String}, pattern::String)

Construct a [`PRS.PatternFreeString`](@ref).

```jldoctest; output = true
using PartialRejectionSampling
PRS.PatternFreeString("ATGTA", ["A", "C", "G", "T"])

# output

PatternFreeString{String}
- pattern = ATGTA
- alphabet = ["A", "C", "G", "T"]
```
"""
function PatternFreeString(pattern::String, alphabet::Vector{String})
    @assert !isempty(pattern)
    @assert !isempty(alphabet)
    if !issubset(string.(unique(pattern)), alphabet)
        throw(DomainError(pattern, "pattern is not fully made of characters from alphabet $(alphabet)"))
    end
    return PatternFreeString{String}(pattern, alphabet)
end

"""
    generate_sample(
        [rng=Random.AbstractRNG,]
        pp::PatternFreeString{T},
        size::Int
    )::T where {T<:String}

Generate a string uniformly at random among all strings made of characters from `pp.alphabet` with no occurence of `pp.pattern`.

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

Generate a string uniformly at random among all strings made of characters from `pp.alphabet` with no occurence of `pp.pattern`, using a tailored version of Partial Rejection Sampling (PRS).

```@example
using PartialRejectionSampling
pp = PRS.PatternFreeString("ATGTA", ["A", "C", "G", "T"])
PRS.generate_sample_prs(pp, 20)
```

**See also**
- [GiAmWe18](@cite), Sections 2.1 and 4
"""
function generate_sample_prs(
    rng::Random.AbstractRNG,
    pp::PatternFreeString{T},
    size::Int
)::T where {T<:String}
    @assert size > 0
    return _pattern_free_string_prs(rng, pp.pattern, pp.alphabet, size)
end

"""
    _generate_pattern_free_string_prs(
        rng::Random.AbstractRNG,
        pattern::String,
        alphabet::Vector{String},
        size::Int
    )::String

Generate a string uniformly at random among all strings made of characters from `alphabet` with no occurence of the pattern `pattern`, using a tailored version of Partial Rejection Sampling (PRS)

**See also**
- [GiAmWe18](@cite), Sections 2.1 and 4
"""
function _pattern_free_string_prs(
    rng::Random.AbstractRNG,
    pattern::String,
    alphabet::Vector{String},
    size::Int
)::String
    if has_common_prefix_suffix(pattern)
        return _pattern_free_string_prs_general(rng, pattern, alphabet, size)
    else
        return _pattern_free_string_prs_extremal(rng, pattern, alphabet, size)
    end
end

## Extremal PRS

function _pattern_free_string_prs_extremal(rng, pattern, alphabet, size)
    str_vec = rand(rng, alphabet, size)
    while true
        str = join(str_vec)
        bad = findall(pattern, str; overlap=false)
        isempty(bad) && return str
        for range_ in bad
            str_vec[range_] .= rand(rng, alphabet, length(range_))
        end
    end
end

## General PRS

function _pattern_free_string_prs_general(rng, pattern, alphabet, size)
    str_vec = Vector{String}(undef, size)
    rand!(rng, str_vec, alphabet)
    tmp_vec = fill("", size)
    while true
        str = join(str_vec)
        bad = findall_overlap(pattern, str)
        isempty(bad) && return str
        for range_ in bad
            @inbounds tmp_vec[range_] .= str_vec[range_]
        end
        res = empty(bad)
        while !isempty(bad)
            B = popfirst!(bad)
            B̄ = _check_extension!(tmp_vec, str_vec, pattern, B)
            if isequal(B, B̄)
                push!(res, B)
            else
                push!(bad, B̄)
            end
        end
        for range_ in res
            @inbounds str_vec[range_] .= rand(rng, alphabet, length(range_))
        end
    end
end

"""
    _check_extension!(
        tmp_vec::Vector{String},
        str_vec::Vector{String},
        pattern::String,
        window::UnitRange{U}
    )::UnitRange{U} where {U<:Int}

Assuming `join(tmp_vec[window]) == join(str_vec[window]) == pattern`, check whether a reassigment of `""` elements from left or right of `tmp_vec[window]` can make `pattern` occur.
If this is the case, the identified `""` elements of `tmp_vec` are set with the corresponding elements from `str_vec` and the *extended* window where pattern can arise is returned.
Otherwise the original `window` is returned.
"""
function _check_extension!(
    tmp_vec::Vector{String},
    str_vec::Vector{String},
    pattern::String,
    window::UnitRange{U}
)::UnitRange{U} where {U<:Int}
    p = length(pattern)
    i, j = first(window), last(window)
    # left looking
    i₋ = max(firstindex(tmp_vec), i-p+1)
    while i₋ < i
        window_ = i₋:min(i₋+p-1, j)
        range_ = _pattern_can_occur_if_reassignment_at(pattern, tmp_vec, window_)
        if !isempty(range_)
            @inbounds tmp_vec[range_] .= str_vec[range_]
            break
        end
        i₋ += 1
    end
    # right looking
    j₊ = min(lastindex(tmp_vec), j+p-1)
    while j₊ > j
        window_ = j₊:max(j₊-p+1, i)
        range_ = _pattern_can_occur_if_reassignment_at(pattern, tmp_vec, window_)
        if !isempty(range_)
            @inbounds tmp_vec[range_] .= str_vec[range_]
            break
        end
        j₊ -= 1
    end
    return i₋:j₊
end

"""
    _pattern_can_occur_if_reassignment_at(
        pattern::String,
        vec::Vector{String},
        window::UnitRange{U}
    )::UnitRange{U} where {U<:Int}

Find the range of indices of `window` where `""` elements of `vec[window]` can be modified to make the resulting `join(vec[window]) == pattern`.
An empty range is returned otherwise.
"""
function _pattern_can_occur_if_reassignment_at(
    pattern::String,
    vec::Vector{String},
    window::UnitRange{U}
)::UnitRange{U} where {U<:Int}
    empty_range = one(U):zero(U)
    (isempty(window) || length(window) != length(pattern)) && return empty_range
    i = first(window)
    f1 = findnext(isempty, vec, i)
    isnothing(f1) && return empty_range
    j = last(window)
    f2 = findprev(isempty, vec, j)
    if startswith(pattern, join(vec[i:f1])) && endswith(pattern, join(vec[f2:j]))
        return f1:f2
    end
    return empty_range
end
