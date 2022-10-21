using .Utils: Uniform, num_possible_bitstrings
using ..Bits: bit_string, AbstractBitVector1
using ..DictTools: _AbstractDict
using ..TupleArrays

export build_dist

"""
    build_dist(edge_vals::AbstractArray, nbits::Integer)

Build a distribution over observed bit sequences from a distribution
of error probabilities over bit sequences.
"""
build_dist(edge_vals::AbstractVector) = _build_dist(edge_vals, length(edge_vals))
build_dist(edge_vals::AbstractVector{<:MaybeNamedTuple}) = _build_dist(edge_vals, num_possible_bitstrings(edge_vals))
build_dist(edge_vals::_AbstractDict) = _build_dist(edge_vals, num_possible_bitstrings(edge_vals))

# edge_pairs returns pairs whose first member is a bitstring (coded as an integer),
# and second is a probabilities.
# Index 1 in Vector corresponds to bitstring 0, so we must adjust
@inline edge_pairs(x::AbstractVector) = ((i-1, v) for (i, v) in pairs(x))
# In a Dict, the key is already correct.
@inline edge_pairs(x::_AbstractDict) = pairs(x)
# This already iterates over pairs
@inline edge_pairs(x::AbstractVector{<:MaybeNamedTuple}) = x

# This is for dense arrays of error strings. Not for realistic (in partic. sparse) data
#function _build_dist(edge_vals::AbstractVector, nbitstrings::Integer)
function _build_dist(edge_vals, nbitstrings::Integer)
    T = typeof(last(first(edge_pairs(edge_vals))))
    prob_dist = zeros(T, nbitstrings)
    new_prob_dist = zeros(T, nbitstrings)
    # start with no error--- concentrate all prob on (0,...,0)
    prob_dist[1] = one(eltype(prob_dist))
    for (edge, edge_p) in edge_pairs(edge_vals)
        iszero(edge_p) && continue
        OMupdate!(new_prob_dist, prob_dist, (edge, edge_p))
        (prob_dist, new_prob_dist) = (new_prob_dist, prob_dist) # swap buffers
    end
    return prob_dist # last swap makes prob_dist the most up to date
end

# Index of 1 represents bitstring (0,...,0).
# So, convert from index to bitstring, do xor, then convert back to index.
@inline _apply_flips(k, v) = typeof(k)(xor(k-1, v) + 1)

"""
    OMupdate!(new_prob_dist::AbstractVector, prob_dist::AbstractVector, (v, p))

Update the probabilities of read bit sequences with hyperedge `v` and
probability `p`.
"""
function OMupdate!(new_prob_dist::AbstractVector, prob_dist::AbstractVector, (v, p))
    for k in eachindex(prob_dist)
       new_prob_dist[k] = (1 - p) * prob_dist[k]
    end
    for k in eachindex(prob_dist)
        i = _apply_flips(k,v)
        new_prob_dist[i] += prob_dist[k] * p
    end
    return nothing
end

"""
    toy_edge_probabilities([T=Float64], nbits::Integer)

Make a toy probability distribution over error patterns, which are tuples of
bits. We represent these bits via the bits in an Int64. In the toy distribution
all patterns have non-zero weight.
"""
toy_edge_probabilities(nbits::Integer) = toy_edge_probabilities(Float64, nbits)
function toy_edge_probabilities(::Type{T}, nbits::Integer) where T
    num_edges = 2 ^ nbits
    dist = Uniform{T}(1e-5, 1e-1)
    edge_vals = rand(dist, num_edges)
    return edge_vals
end
