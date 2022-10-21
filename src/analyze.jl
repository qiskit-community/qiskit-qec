module Analyze

using ..Bits
using ..Bits._Bits: BitVector1Mask, AbstractBitVector1, bits, bit
using ..Utils: min_bits
using Dictionaries
using ..DictTools: count_map, add_ncounts!, hist_to_dist

export bit_touch_count, all_bits_touched, error_string_info,
    sampled_error_info,
    freq_one_count, prob_one_error, prob_two_errors, prob_n_errors

export matchat, patcountmap

_torange(n::Integer, m::Integer) = n:(n+m-1)

"""
    matchat(str::AbstractString, start, x::Bits._Bits.BitVector1Mask)

Return `true` if the bitstring `str` is a "substring" of `x` at index `start`.
That is, bits `start : (start+length(str)-1)` are equal to `str`.
"""
function matchat(str::AbstractString, start, x::Bits._Bits.BitVector1Mask)
    mnum = bits(str)
    return mnum == x[_torange(start, length(mnum))]
end

"""
    matchat(tests::Union{AbstractArray, Tuple}, x)

Return `true` if each test in `test` returns `true`.

Elements of `test` are `Tuples` of the form `(str, start)` which are passed

"""
function matchat(tests::Union{AbstractArray,Tuple}, x::Bits._Bits.BitVector1Mask)
    return all(y -> matchat(y..., x) for y in tests)
end

# This would be more efficient if we precomputed the match strings
function matchat(tests::Union{AbstractArray,Tuple})
    return x -> all(matchat(y..., x) for y in tests)
end

"""
    matchat(str::AbstractString, start)

Like `matchat(str, start, x)`, but return a function that is then
applied to `x`. In other words, curry the first two arguments.

# Examples
```julia-repl
julia> eprobs = read_hyper_edges(NamedTupleVector, preprocessed_filename);
julia> filter(x -> matchat("01", 1)(x[1]), eprobs)
```
The following is ten times faster (prbly because our indexing code is inefficient).
But, this only returns the bitstrings, not the probs.
But, the example above dispatches to the `filter(f, a::AbstractArray)` method,
which is slow.
```julia-repl
julia> filter(matchat("01", 1), eprobs(:bitstr))
```
Yes, the following is slow-ish
```julia-repl
julia> findall(x -> matchat("01", 1)(x[1]), eprobs)
```
"""
function matchat(str::AbstractString, start)
    mnum = bits(str)
    rng = _torange(start, length(mnum))
    return x -> mnum == x[rng]
end

"""
    patcountmap(r::UnitRange, estrings)

Make a countmap of bit patterns in the collection of bitstrings `estrings`
where testing for equality is restricted to the bit-indices given by `r`.
That is `<010111>` and `<000110>` are counted in the same bin for `r=2:3`,
because their second and third bits are the same.

# Examples
```julia-repl
julia> patcountmap(1:3, eprobs(:bitstr))
5-element Dictionary{BitVector1Mask{UInt128}, Int64}
 <000> │ 196
 <001> │ 48
 <100> │ 30
 <010> │ 8
 <110> │ 6
```
"""
patcountmap(r::UnitRange, estrings) = sort(count_map(x[r] for x in estrings); rev=true)
patcountmap(n, m, estrings) = patcountmap(_torange(n, m), estrings)

"""
    bittouch(r::UnitRange, x)

Return `true` if any bits in the range `r` are equal to `1`.
That is, does this error string "touch" or flip bits.
"""
bittouch(r::UnitRange, x) = count_ones(x[r]) > 0

"""
    bittouch(r::UnitRange)

Curried version of `bittouch(r::UnitRange, x)`.
"""
bittouch(r::UnitRange) = x -> bittouch(r, x)
bittouch(n::Integer, m::Integer, args...) = bittouch(_torange(n, m), args...)
bittouch(n::Integer, m::Integer) = bittouch(_torange(n, m))

###
### Analyze error probability map (or just the error strings).
### Not the observed bitstrings.
###

"""
    bit_touch_count(bstrings)

For each bit position `i`, count the number elements of `bstrings` that have `1` that
position.
`bstrings` are binary strings represented as bitstype integers.
A common prefix of leading zeros is ignored.
"""
function bit_touch_count(bstrings)
    nbits = maximum(min_bits, bstrings)
    return [count(isone(bit(bstrings[j], i)) for j in eachindex(bstrings)) for i in 1:nbits]
end

"""
    all_bits_touched(bstrings)

Return `true` if each bit position is `1` for at least one element of `bstrings`.
Equivalently, return `false` if for some bit position `i`, none of `bstrings`
has that bit set.
A common prefix of leading zeros is ignored.
"""
all_bits_touched(bstrings) = all(!iszero, bit_touch_count(bstrings))


"""
    error_string_info(error_probs)

Compute some information on the dict of error probabalities `error_probs`.

pr_all_occur -- (prob that all errors occur, decimal exponent of prob)
pr_n_occur -- tuple of probs that exactly (0, 1, 2, 3, 4) errors occur.
all_bits_touched -- true if each bit is flipped by at least one error.
                    If false, then this bit is irrelevant.
                    See `bit_flip_count()`.
"""
function error_string_info(error_probs)
    bstrings = keys(error_probs)
    probs = values(error_probs)
    pr_all_occur =
        setprecision(2^13) do
            p_all = prod(big.(probs))
            _all = (Float64(p_all), floor(Int, log10(p_all)))
            return _all
        end

    pr_n_occur = prob_n_errors(error_probs)

    _all_bits_touched = all_bits_touched(bstrings)

    return (; pr_all_occur, pr_n_occur, all_bits_touched=_all_bits_touched)
end

"""
    prob_n_errors(error_probs, n=4)

Return Tuple of probabilities that exactly 0, 1, 2, 3, 4, and 5,
errors occur.

`error_probs` is a map of errors (bitstrings) to probabilities. Each error
occurs independently with probability given by the map.
Passing `n` satisfying `n in 0:5` limits the number of probabilities returned.
Increasing `n` by one increases run time by a factor of about `length(error_probs)`.
The factor may be, for instance, 200.
For n>4, it may be sufficient to get results to a few digits of precision using
`one_sample_num_errors!`
"""
function prob_n_errors(error_probs, n=4)
    n in 0:5 || throw(DomainError("n must be in (0,1,2,3,4,5), got $n"))
    probs = (collect(values(error_probs))...,) # using Tuple is faster below
    m = length(probs)
    fac = i -> probs[i] / (1 - probs[i])

    p_none = prod(1 - x for x in probs)
    f0 = () -> p_none
    f1 = () -> p_none * sum(fac(i) for i in 1:m)
    f2 = () ->
        (indpairs = ((i, j) for i in 1:m for j in (i+1):m);
        p_none * sum(fac(i) * fac(j) for (i, j) in indpairs))

    f3 = () ->
        (indtrips = ((i, j, k) for i in 1:m for j in (i+1):m for k in (j+1):m);
        p_none * sum(fac(i) * fac(j) * fac(k) for (i, j, k) in indtrips))

    f4 = () ->
        (indquad = ((i, j, k, l) for i in 1:m for j in (i+1):m for k in (j+1):m for l in (k+1):m);
        p_none * sum(fac(i) * fac(j) * fac(k) * fac(l) for (i, j, k, l) in indquad))

    f5 = () ->
        (ind5 = ((i, j, k, l, m) for i in 1:m for j in (i+1):m for k in (j+1):m for l in (k+1):m for m in l:m);
        p_none * sum(fac(i) * fac(j) * fac(k) * fac(l) * fac(m) for (i, j, k, l, m) in ind5))

    pfuncs = (f0, f1, f2, f3, f4, f5)
    return ((pfuncs[i]() for i in 1:(n+1))...,)
end

###
### Analyze countmap of sampled outcomes
###

"""
    sampled_error_info(samp_map)

Info on `samp_map`, whose keys are bitstrings (as integers) and values
are the number of times the string was observed. `samp_map` is generated
via Monte Carlo. Each sample is generated by randomly applying each error with its
specified probability. Note that several error combinations may give rise to
the same observed bit string. This information is lost in `samp_map`.
"""
function sampled_error_info(samp_map)
    num_counts = sum(values(samp_map))
    density = length(samp_map) / num_counts
    freq_one_count_no_weights = hist_to_dist(freq_one_count(samp_map; weighted=false))
    freq_one_count_weights = hist_to_dist(freq_one_count(samp_map; weighted=true))
    return (; freq_one_count_no_weights, freq_one_count_weights, num_counts, density)
end

function freq_one_count(samp_map; weighted=false)
    if weighted
        _map = Dictionary{Int,Int}()
        for (bstring, _count) in pairs(samp_map)
            add_ncounts!(_map, count_ones(bstring), _count)
        end
        return sortkeys!(_map)
    else
        return sortkeys!(count_map(Dictionary, count_ones.(keys(samp_map))))
    end
end



end # module Analyze
