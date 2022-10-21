module Utils

import Printf
import Dictionaries
using ..DictTools: _AbstractDict
using ..TupleArrays
using ..HeavyHex: MaybeNamedTuple
using ..Bits: min_bits

# TODO: export fewer things when we aren't testing so much
export Uniform, num_possible_bitstrings, to_dense_vector

import Random

# Loading Distributions still takes 0.5s, so we roll our own distribution...
# We could use Float64 instead of T for params. Maybe computation is more accurate that way.
"""
    Uniform{[T<:AbstractFloat = Float64]}(min, max)

Uniform distribution over (`min`, `max`). A random sampler is defined so that
`rand(Uniform(.1, .2), dims)` works and is efficient.

# Examples
```julia-repl
julia> rand(Uniform(.1, .2))
0.135336264662535

julia> rand(Uniform(Float32(.1), Float32(.2)))
0.15677246f0
```
"""
struct Uniform{T<:AbstractFloat}
    min::T
    max::T
end
Uniform(args...) = Uniform{Float64}(args...) # default is Float64
Base.eltype(::Type{Uniform{T}}) where T = T

# Define the sampler, so that methods of `rand` work.
function Random.rand(rng::Random.AbstractRNG, s::Random.SamplerTrivial{Uniform{T}}) where T
    return s[].min + (s[].max-s[].min)*rand(rng, T)
end

# Probably don't need this ? Too simple?
# """
#     min_bits(d::_AbstractDict)

# Find the longest bitstring in the binary representation of the
# integer valued keys of `d`, discounting leading zeros.
# """
# min_bits(d::_AbstractDict) = maximum(min_bits, keys(d))

# """
#     min_bits(d::AbstractVector)

# PROBABLY DONT WANT TO USE THIS, use maximum(min_bits, v) instead?
# Find the longest bitstring in the binary representation of the
# integer valued keys of `d`, discounting leading zeros.
# """
# min_bits(v::AbstractVector) = maximum(min_bits, v)

"""
    num_possible_bitstrings(d::_AbstractDict)

Return the number of possible bitstrings of length equal to the longest bitstring in the keys of `d`.
Leading zeros in the keys of `d` are ignored.

For example if the longest bitstring in `keys(d)` has length `n`, then return `2^n`.
"""
function num_possible_bitstrings(d) #::_AbstractDict)
    maxbits = 30
    nbits = maximum(min_bits, _get_bitstrs(d))
    nbits > maxbits &&
        error("Max length of bitstring is $maxbits. Won't try to allocate 8 * 2^$nbits = $(8 * 2^nbits) bytes = $(_gbfrombits(nbits)) GB.")
    return 2^nbits
end

_get_bitstrs(d::_AbstractDict) = keys(d)
_get_bitstrs(d::NamedTupleArray) = d(1)
# I would like to restrict this more, but thats a pain
_get_bitstrs(d::AbstractVector{<:MaybeNamedTuple}) = (first(x) for x in d)

_gbfrombits(n) = _gbstring(8 * 2.0^n)

_gbstring(x) = Printf.@sprintf("%.2f", x / 1e9)

"""
    to_dense_vector(dict::_AbstractDict{<:Integer, ValT}) where {ValT}

Convert `dict` to `v::Vector{ValT}`, assuming the keys of `dict` are integers
representing the indices of `v`. The keys are incremented by one. For
example, a key of `0` is the first index in `v`.
"""
function to_dense_vector(dict::_AbstractDict{<:Integer, ValT}) where {ValT}
    arr = zeros(ValT, num_possible_bitstrings(dict))
    for (ind, val) in dict
        arr[ind + 1] = val
    end
    return arr
end

end # module Utils
