module SmallDicts

# I think this is kind of obsolete
# Mostly superceded by TupleArrays; a few features are not implemented there yet.

using Base: Ordering, Forward

export SmallDict, SmallDictNamed, small_dict, keyname, valuename, tuples, NamedTupleVector, checklength

function small_dict end

# This does not work
if @isdefined AbstractDictionary
    SmallDict(dict::AbstractDictionary, sorted=false, presorted=false) = dict_to_smalldict(dict, sorted, presorted)
    Dictionary(sd::SmallDict) = Dictionary(keys(sd), values(sd))
    small_dict(d::AbstractDictionary; keyname=nothing, valuename=nothing) =
        small_dict(keys(d), values(d); keyname=keyname, valuename=valuename)
end

abstract type AbstractSmallDict{K,V,S,N} <:  AbstractDict{K,V} end

#const _SmallDict = SmallDict
const _SmallDict = AbstractSmallDict

"""
    SmallDict{...}

Immutable struct that behaves like a `Dict` when iterating over `pairs(::SmallDict)`.
Iterating is 100 times faster than over `Dict`. Data is stored as two `NTuples{N}`.
`N` should not be too large for good performance. How big ? 100,200,300?

Converting a small dictionary to `SmallDict` will allow you to iterate over
key-value pairs much more efficiently. Sometimes we iterate repeatedly over
a dictionary doing very little work at each iteration. When generating Monte Carlo
data, this may yield a large savings in time.
"""
struct SmallDict{K,V,S,N} <: AbstractSmallDict{K,V,S,N}
    _keys::NTuple{N,K}
    vals::NTuple{N,V}
end

struct SmallDictNamed{K,V,S,N,KName,VName} <: AbstractSmallDict{K,V,S,N}
    _keys::NTuple{N,K}
    vals::NTuple{N,V}
end

_tsort(args...) = sort(args...)
_tsort(t::Tuple) = (sort!([t...])...,)

function SmallDictNamed(kname::Symbol, vname::Symbol, _inkeys, _invals, sorted::Bool=false, presorted=false)
    n = length(_inkeys)
    length(_invals) == n || throw(DimensionMismatch("Keys and values have different lengths."))
    inkeys = sorted && ! presorted ? _tsort(_inkeys) : _inkeys
    invals = sorted && ! presorted ? _tsort(_invals) : _invals
    _keys = (inkeys...,)
    vals = (invals...,)
    return SmallDictNamed{eltype(_keys), eltype(vals), sorted, n, kname, vname}(_keys, vals)
end

function Base.getproperty(sd::SmallDictNamed, prop::Symbol)
    prop === keyname(sd) && return keys(sd)
    prop === valuename(sd) && return values(sd)
    return getfield(sd, prop)
end

Base.propertynames(sd::SmallDictNamed) = (keyname(sd), valuename(sd))

function small_dict(args...; keyname=nothing, valuename=nothing)
    (knothing, vnothing) = (isnothing(keyname), isnothing(valuename))
    knothing == vnothing || throw(ErrorException("Both or neither name for keys and vals must be supplied"))
    isnothing(keyname) &&  return SmallDict(args...)
    return SmallDictNamed(keyname, valuename, args...)
end

"""
    SmallDict{N,TK,TV}(dict::_AbstractDict)

Convert a smallish `dict` to an (immutable) `SmallDict`.
"""
function SmallDict(_inkeys, _invals, sorted::Bool=false, presorted=false)
    n = length(_inkeys)
    length(_invals) == n || throw(DimensionMismatch("Keys and values have different lengths."))
    inkeys = sorted && ! presorted ? _tsort(_inkeys) : _inkeys
    invals = sorted && ! presorted ? _tsort(_invals) : _invals
    _keys = (inkeys...,)
    vals = (invals...,)
    return SmallDict{eltype(_keys), eltype(vals), sorted, n}(_keys, vals)
end

dict_to_smalldict(dict, sorted=false, presorted=false) = SmallDict(keys(dict), values(dict), sorted, presorted)

SmallDict(dict::AbstractDict, sorted=false, presorted=false) = dict_to_smalldict(dict, sorted, presorted)
Dict(sd::_SmallDict) = Dict(zip(keys(sd), values(sd)))

keyname(d::SmallDictNamed{<:Any,<:Any,<:Any,<:Any,KName,<:Any}) where KName = KName
valuename(d::SmallDictNamed{<:Any,<:Any,<:Any,<:Any,<:Any,VName}) where VName = VName

Base.keys(sd::_SmallDict) = sd._keys
Base.values(sd::_SmallDict) = sd.vals
@inline Base.pairs(sd::_SmallDict{<:Any,<:Any,S, N}) where {S,N} = sd
@inline Base.length(::_SmallDict{<:Any,<:Any,S, N}) where {S,N} = N
@inline Base.getindex(sd::_SmallDict, inds...) = (keys(sd)[inds...], values(sd)[inds...])
@inline Base.lastindex(sd::_SmallDict{<:Any,<:Any,S,N}) where {S,N} = keys(sd)[N]
@inline Base.firstindex(sd::_SmallDict) = keys(sd)[1]
_pair(sd::_SmallDict, i::Integer) = Pair(keys(sd)[i], values(sd)[i])
@inline Base.iterate(sd::_SmallDict, i=1) = i > length(sd) ? nothing : (_pair(sd, i), i + 1)

# propagate probably does nothing
@inline Base.@propagate_inbounds function _pair(sd::SmallDictNamed, i::Integer)
    tup = (keys(sd)[i], values(sd)[i])
    return NamedTuple{(keyname(sd), valuename(sd)), typeof(tup)}(tup)
end

@inline function getposition(tup::Union{Tuple, AbstractVector}, _key)::Int
    i::Int = 1
    for n in tup
        _key == n && return i
        i += 1
    end
    return i
end

@inline getposition(sd::_SmallDict{<:Any, <:Any, false}, _key)::Int = getposition(keys(sd), _key)
@inline getposition(sd::_SmallDict{<:Any, <:Any, true}, _key)::Int = searchsortedfirst(keys(sd), _key)

# This is faster if the index is high enough
@inline function Base.getindex(sd::_SmallDict, _key)
    ind = getposition(sd, _key)
    ind > length(sd) && throw(KeyError(_key))
    return values(sd)[ind]
end

# index of the first value of vector a that is greater than or equal to x;
# returns length(v)+1 if x is greater than all values in v.
function Base.searchsortedfirst(v::Tuple, x, lo::T, hi::T, o::Ordering) where T<:Integer
    u = T(1)
    lo = lo - u
    hi = hi + u

    @inbounds while lo < hi - u
        m = Base.Sort.midpoint(lo, hi)
        if Base.Sort.lt(o, v[m], x)
            lo = m
        else
            hi = m
        end
    end
    return hi
end

Base.searchsortedfirst(v::Tuple, x;
           lt=isless, by=identity, rev::Union{Bool,Nothing}=nothing, order::Ordering=Forward) =
            Base.searchsortedfirst(v,x,Base.ord(lt,by,rev,order))

Base.searchsortedfirst(v::Tuple, x, o::Ordering) = (inds = Base.OneTo(length(v)); searchsortedfirst(v,x,first(inds),last(inds),o))



end # module SmallDicts
