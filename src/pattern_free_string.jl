"""
    PatternFreeString{T<:String} <: AbstractPointProcess{T}

Container with fields `alphabet` and `pattern` used to generate a string made of characters from `alphabet` avoiding the prescribed `pattern`
"""
struct PatternFreeString{T<:String} <: AbstractPointProcess{T}
    alphabet::Vector{T}
    pattern::T
end

"""
    PatternFreeString(alphabet::Vector{String}, pattern::String)

Construct a [`PatternFreeString`](@ref)
"""
function PatternFreeString(alphabet::Vector{String}, pattern::String)
    @assert !isempty(alphabet)
    @assert !isempty(pattern)
    return PatternFreeString{String}(alphabet, pattern)
end

"""
    generate_sample(
        psf::PatternFreeString{T},
        size::Int;
        rng=-1
    )::T where {T<:String}

Generate a string uniformly at random among all strings made of characters from `psf.alphabet` with no occurence of the pattern `psf.pattern`.

Default sampler is [`generate_sample_prs`](@ref)
"""
function generate_sample(
        psf::PatternFreeString{T},
        size::Int;
        rng=-1
)::T where {T<:String}
    return generate_sample_prs(psf, size; rng=rng)
end

"""
    generate_sample_prs(
        psf::PatternFreeString{T},
        size::Int;
        rng=-1
    )::T where {T<:String}

Generate a string uniformly at random among all strings made of characters from `psf.alphabet` with no occurence of the pattern `psf.pattern`, using a tailored version of Partial Rejection Sampling (PRS) derived by [GiAmWe18](@cite)
"""
function generate_sample_prs(
        pfs::PatternFreeString{T},
        size::Int;
        rng=-1
)::T where {T<:String}
    @assert size > 0
    return _generate_sample_pattern_free_string_prs(pfs.alphabet, pfs.pattern, size; rng=rng)
end

@doc raw"""
    generate_pattern_free_string_prs(
        alphabet::Vector{T},
        pattern::T,
        size::Int;
        rng=-1
    )::T where {T<:AbstractString}

Generate a string uniformly at random among all strings made of characters from `alphabet` with no occurence of the pattern `pattern`, using a tailored version of Partial Rejection Sampling (PRS) derived by [GiAmWe18](@cite)
"""
function _generate_sample_pattern_free_string_prs(
    alphabet::Vector{T},
    pattern::T,
    size::Int;
    rng=-1
)::T where {T<:AbstractString}

    @assert size > 0
    @assert !isempty(alphabet)
    @assert !isempty(pattern)
    rng = getRNG(rng)

    pref_suff = find_prefix_suffix(pattern)

    pfs = fill("", size)
    resample_indices = Set(1:size)

    while !isempty(resample_indices)
        generate_sample!(pfs, resample_indices, alphabet; rng=rng)
        resample_indices = find_characters_to_resample(pfs, pattern, pref_suff)
    end
    return join(pfs)
end

find_prefix_suffix(s::String) = [i for i in 1:div(length(s), 2) if s[1:i] == s[end-i+1:end]]

"""
    generate_sample!(
        string_vec::Vector{T},
        indices,
        alphabet::Vector{T};
        rng=-1
    ) where {T<:AbstractString}

Generate a character uniformly at random from `alphabet` at positions prescribed by `indices` in `string_vec`
"""
function generate_sample!(
    string_vec::Vector{T},
    indices,
    alphabet::Vector{T};
    rng=-1
) where {T<:AbstractString}
    rng = getRNG(rng)
    for i in indices
        string_vec[i] = rand(rng, alphabet)
    end
end

"""
    find_bad_ranges(
        pattern::T,
        string::T
    )::Vector{UnitRange} where {T<:AbstractString}

Identify where `pattern` occur in `sting` and return the corresponding ranges of indices.
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
Return the indices of the variables (here characters) involved in the corresponding events.

This function is used as a subroutine of the grid PRS methodology of [MoKr20](@cite), see [`generate_sample_grid_prs`](@ref)

**See also**

[`find_bad_ranges`](@ref)
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
