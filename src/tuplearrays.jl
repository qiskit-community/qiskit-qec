module TupleArrays

export AbstractTupleArray, TupleArray, TupleVector,
    AbstractNamedTupleArray, NamedTupleArray, NamedTupleVector
export TupleMatrix, TupleVector, AbstractTupleVector, AbstractTupleMatrix,
    NamedTupleVector, NamedTupleMatrix
export tupgetdata, tupgetindex, tuplength, getnames, tuptype, getcol, tupget_unnamed_data

# temporary export
export _striptype, ensurevectorcols, hastuplecol, hasvectorcol, alltuplecol, ensuretuplecols

###
### AbstractTupleArray
###

"""
    AbstractTupleArray{Ttup, Nd, Nt} <: AbstractArray{Ttup, Nd}

An abstract array of `Tuple`s. In the concrete subtypes, the
data is not stored as `Tuple`s, but rather as `Nt` `Arrays`, where
`Nt` is the length of each `Tuple`.

`Ttup` is the data type of the array. It is `<:Tuple`.
`Nd` is the dimension of the array.
`Nt` is the length of each `Tuple` in the array.

# Examples

If `Nt=4` and `Nd=2`, then the array is a matrix of `4-Tuples`.
"""
abstract type AbstractTupleArray{Ttup, Nd, Nt} <: AbstractArray{Ttup, Nd} end
const AbstractTupleVector{Ttup} = AbstractTupleArray{Ttup, 1, Nt} where {Ttup, Nt}
const AbstractTupleMatrix{Ttup} = AbstractTupleArray{Ttup, 2, Nt} where {Ttup, Nt}

"""
    tuplength(ta::AbstractTupleArray)

Length  of each `Tuple` contained in `ta`.
"""
tuplength(ta::AbstractTupleArray{<:Any,<:Any,Nt}) where Nt = Nt

"""
    tuptype(ta::AbstractTupleArray)

Type of the `Tuple`s contained in `ta`.
"""
tuptype(ta::AbstractTupleArray{Ttup}) where Ttup = Ttup

tupgetdata(ta::AbstractTupleArray) = getfield(ta, :_data)

tupget_unnamed_data(ta::AbstractTupleArray) = tupgetdata(ta)

Base.ndims(ta::AbstractTupleArray{<:Any, Nd}) where Nd = Nd

"""
    hastuplecol(ta::AbstractTupleArray)

Return `true` if any column in `ta` is a `Tuple`.

The elements of `ta` are always `Tuple`s. But, the columns may
other containers.
"""
hastuplecol(ta::AbstractTupleArray) = any(x->isa(x,Tuple), tupgetdata(ta))

#
alltuplecol(ta::AbstractTupleArray) = all(x->isa(x,Tuple), tupgetdata(ta))
hasvectorcol(ta::AbstractTupleArray) = any(x->isa(x,AbstractVector), tupgetdata(ta))

"""
    getcol(ta::AbstractTupleArray, i::Integer)

Return the `i`th column of `ta`, that is the `i`th element of each
`Tuple` in `ta`.

# Examples
```julia-repl
julia> v = NamedTupleVector((:a, :b), (3, 2, 1), (6, 5, 4))
3-element NamedTupleVector{Tuple{Tuple{Int64, Int64, Int64}, Tuple{Int64, Int64, Int64}}, 2, (:a, :b)}:
 (a = 3, b = 6)
 (a = 2, b = 5)
 (a = 1, b = 4)

julia> getcol(v, 1)
(3, 2, 1)

julia> getcol(v, :b)
(6, 5, 4)
```
"""
getcol(ta::AbstractTupleArray, i) = tupgetdata(ta)[i]

# Make AbstractTupleArray callable. Use for indexing into data
(ta::AbstractTupleArray)(i) = getcol(ta, i)
(ta::AbstractTupleArray)(i, inds...) = getcol(ta, i)[inds...]

Base.getproperty(ta::AbstractTupleArray, s::Symbol) = getcol(ta, s)


"""
    ensurevectorcols(ta)

Copy the tuple of containers `fieldtup` converting any `Tuple` to a `Vector`.
This is probably not robust. Hard to say what needs to be converted.

One reason you might use `ensurevectorcols` is to mutate the result.
"""
ensurevectorcols(fieldtup) = ((isa(x,Tuple) ? collect(x) : x for x in fieldtup)...,)
ensurevectorcols(ta::AbstractTupleArray) = !hastuplecol(ta) ? ta :
    _striptype(typeof(ta))(ensurevectorcols(tupgetdata(ta))...,)

"""
    ensuretuplecols(ta::AbstractTupleArray)

If any columns in `ta` are not `Tuple`s return a copy with all columns
converted to `Tuple`s. Otherwise return `ta`.

If `ta` is not too large, then iterating over and indexing into the returned
array can be much faster.
"""
ensuretuplecols

# TODO: what is efficient ? Tuple(x), (x...) ?
ensuretuplecols(fieldtup) = ((isa(x,Tuple) ? x : Tuple(x) for x in fieldtup)...,)
ensuretuplecols(ta::AbstractTupleArray) = alltuplecol(ta) ? ta :
    _striptype(typeof(ta))(ensuretuplecols(tupgetdata(ta))...,)


# Worked around this by converting AbstractTupleArray with tuple data
# to one with vector data when calling `show`. So we don't need the
# piracy below.
# piracy here. There is some bug in printing that makes
# this necessary. Doing this:
# show(io, MIME"text/plain"(), TupleArray((1,2),(3,4))
# attempts getindex with two indices instead of one. Don't know why or how
Base.getindex(tup::Tuple, i::Integer, j::Integer) = getindex(tup, i)

# The following kinda works, and makes the pirate-hack above uneccessary.
# But, the below will print the wrong type if it converts Tuples to Vectors.
# So we use the pirate hack, instead.
# function Base.show(io::IO, ::MIME"text/plain", ta::AbstractTupleArray)
#     tav = hastuplecol(ta) ? ensurevectorcols(ta) : ta
#     nd = ndims(ta)
#     invoke(show, Tuple{IO,MIME"text/plain", AbstractArray{<:Any,nd}}, io, MIME"text/plain"(), tav)
# end

_tupgetindex(ta::AbstractTupleArray, i...) = ((x[i...] for x in tupgetdata(ta))...,)
tupgetindex(args...) = _tupgetindex(args...)

_size(args...) = Base.size(args...)
_size(t::Tuple) = (length(t),) # Preclude using Tuples for higher dim arrays
# These have length, but you can't index into them, so no good.
# _size(t::Base.ValueIterator) = (length(t),)
# _size(t::Base.KeySet) = (length(t),)

function Base.size(ta::AbstractTupleArray)
    data = tupgetdata(ta)
    isempty(data) && return ((0 for _ in 1:ndims(ta))...,)
    return _size(first(data))
end

Base.length(ta::AbstractTupleArray) = prod(size(ta))

Base.getindex(ta::AbstractTupleArray, i...) = tupgetindex(ta, i...)
Base.getindex(ta::AbstractTupleArray, ind::AbstractVector) = typeof(ta)(tupgetindex(ta, ind)...)

function Base.setindex!(ta::AbstractTupleArray{<:Ttup}, val::Ttup, i::Integer...) where {Ttup}
    data = tupgetdata(ta)
    for j in 1:length(data)
        @inbounds dj = data[j] # This does not help
        @inbounds vj = val[j]
        setindex!(dj, vj, i...)
    end
    return val
end

# TODO: Fix method to check that _Dt is correct if passed
function _construct(arrays, dims=-1, _Nt=-1, _Dt=nothing)
    isempty(arrays) && throw(MethodError(TupleArray, ()))
    asize = _size(first(arrays))
    all(x -> _size(x) == asize, arrays) || throw(DimensionMismatch("Elements of Tuples must have the same dimension"))
    Nd = length(asize)
    dims >= 0 && dims != Nd &&
        throw(DimensionMismatch("Can't construct a $(dims)-dimensional TupleArray with $(Nd)-dimensional data"))
    Dt = typeof(arrays)
    Ttup = Tuple{eltype.(arrays)...}
    Nd = length(asize)
    Nt = length(arrays)
    _Nt >=0 && Nt != _Nt &&
        throw(DimensionMismatch("Can't construct a Tuples of length $_Nt with $Nt input arrays."))
    return (Ttup, Nd, Nt, Dt)
end

struct TupleArray{Ttup, Nd, Nt, Dt} <: AbstractTupleArray{Ttup, Nd, Nt}
    _data::Dt

    TupleArray(arrays...) = new{_construct(arrays)...}(arrays)
    TupleArray{Ttup}(arrays...) where {Ttup} = new{_construct(arrays)...}(arrays)
    TupleArray{Ttup,Nd}(arrays...) where {Ttup, Nd} = new{_construct(arrays, Nd)...}(arrays)
    TupleArray{Ttup,Nd,Nt}(arrays...) where {Ttup, Nd, Nt} = new{_construct(arrays, Nd, Nt)...}(arrays)
    TupleArray{Ttup,Nd,Nt,Dt}(arrays...) where {Ttup, Nd, Nt, Dt} = new{_construct(arrays, Nd, Nt)...}(arrays)
end

TupleArray(d::AbstractDict) = TupleArray(collect(keys(d)), collect(values(d)))

const TupleVector{Ttup, Nt} = TupleArray{Ttup, 1, Nt}
const TupleMatrix{Ttup, Nt} = TupleArray{Ttup, 2, Nt}

####
#### Convert
####

# I wrote an Array method, but it was slower than the fallback method for collect
Base.Array(ta::AbstractTupleArray) = error("Use collect(ta) instead of Array(ta)")
Base.Vector(ta::AbstractTupleArray) = error("Use collect(ta) instead of Vector(ta)")

function Base.copy!(dst::T, src::T) where T <: AbstractTupleArray
    for i in 1:tuplength(src)
        copy!(getcol(dst, i), getcol(src, i))
    end
    return dst
end

for f in (:copy, :similar, :empty)
    @eval (Base.$f)(ta::AbstractTupleArray) = typeof(ta)((($f)(getcol(ta, i)) for i in 1:tuplength(ta))...)
end

function Base.sort!(ta::AbstractTupleArray; col=1, kwargs...)
    ix = sortperm(getcol(ta, col); kwargs...)
    for i in 1:tuplength(ta)
        permute!(getcol(ta, i), ix)
    end
    return ta
end

Base.sort(ta::AbstractTupleArray; kwargs...) = sort!(copy(ta); kwargs...)

###
### AbstractNamedTupleArray
###

abstract type AbstractNamedTupleArray{Ttup, Nd, Nt, Tnames} <: AbstractTupleArray{Ttup, Nd, Nt} end

# Allows a common interface
tupget_unnamed_data(ta::AbstractNamedTupleArray) = values(tupgetdata(ta))

# TODO: should we define keys ?
getnames(::AbstractNamedTupleArray{<:Any,<:Any,<:Any,Tnames}) where Tnames = Tnames

Base.propertynames(ta::AbstractNamedTupleArray) = getnames(ta)

Base.collect(ta::AbstractNamedTupleArray) = [ta...]

function ensurevectorcols(ta::AbstractNamedTupleArray)
    return _striptype(typeof(ta))(getnames(ta), ensurevectorcols(tupgetdata(ta))...,)
end

function ensuretuplecols(ta::AbstractNamedTupleArray)
    return _striptype(typeof(ta))(getnames(ta), ensuretuplecols(tupgetdata(ta))...,)
end

function tupgetindex(ta::AbstractNamedTupleArray, i...)
    tup = _tupgetindex(ta, i...)
    names = getnames(ta)
    return NamedTuple{names, typeof(tup)}(tup)
end

# May be a clever way to avoid defining this. But, it is ok
Base.getindex(ta::AbstractNamedTupleArray, ind::AbstractVector) = NamedTupleArray(tupgetindex(ta, ind))

function _check_name_length(names, arrays)
    length(names) == length(arrays) ||
        throw(DimensionMismatch(
            "Number of names $(length(names)) differ from number of columns $(length(arrays))."))
    return true
end



# TODO: Many or most of these do not need to be internal constructors
struct NamedTupleArray{Ttup, Nd, Nt, Dt, Tnames} <: AbstractNamedTupleArray{Ttup, Nd, Nt, Tnames}
    _data::NamedTuple{Tnames, Dt}

    function NamedTupleArray(ntup::NamedTuple)
        return NamedTupleArray(keys(ntup), ntup...)
    end

    function NamedTupleArray(names::NTuple, arrays...)
        _check_name_length(names, arrays)
        return new{_construct(arrays)..., names}(NamedTuple{names, typeof(arrays)}(arrays))
    end

    function NamedTupleArray{Ttup,Nd}(names::NTuple, arrays...) where {Ttup, Nd}
        _check_name_length(names, arrays)
        return new{_construct(arrays, Nd)..., names}(NamedTuple{names, typeof(arrays)}(arrays))
    end

    function NamedTupleArray{Ttup,Nd,Nt}(names::NTuple, arrays...) where {Ttup, Nd, Nt}
        _check_name_length(names, arrays)
        return new{_construct(arrays, Nd, Nt)..., names}(NamedTuple{names, typeof(arrays)}(arrays))
    end
end

#function NamedTupleArray{<:Any, <:Any, <:Any, <:Any, Tnames}(data...) where Tnames
function NamedTupleArray{Ttup, Nd, Nt, Dt, Tnames}(data...) where  {Ttup, Nd, Nt, Dt, Tnames}
    return NamedTupleArray{Ttup, Nd, Nt}(Tnames, data...)
end

# Allows NamedTupleArray(a = [1,2], b = [3,4]) , etc.
NamedTupleArray(; kwargs...) = NamedTupleArray(keys(kwargs), values(kwargs)...)

NamedTupleArray(names::NTuple, d::AbstractDict) = NamedTupleArray(names, collect(keys(d)), collect(values(d)))

"""
    NamedTupleArray(names::NTuple, ta::AbstractTupleArray)

Construct a NamedTupleArray from `ta` using `names`.
"""
NamedTupleArray(names::NTuple, ta::AbstractTupleArray) = NamedTupleArray(names, tupgetdata(ta)...)

"""
    TupleArray(nta::NamedTupleArray)

Convert `nta` to a `TupleArray` (with no names).
"""
TupleArray(nta::NamedTupleArray) = TupleArray(tupgetdata(nta)...)

const NamedTupleVector{Ttup, Nt} = NamedTupleArray{Ttup, 1, Nt}
const NamedTupleMatrix{Ttup, Nt} = NamedTupleArray{Ttup, 2, Nt}

TupleVector(nta::NamedTupleArray) = TupleArray(nta)
TupleVector(args...) = TupleArray{typeof(args),1}(args...)
TupleMatrix(args...) = TupleArray{typeof(args),2}(args...)

NamedTupleVector(names, args...) = NamedTupleArray{typeof(args),1}(names, args...)
NamedTupleVector(nt::NamedTuple) = NamedTupleArray(nt) # TODO verify dimension

NamedTupleMatrix(names, args...) = NamedTupleArray{typeof(args),2}(names, args...)

for _type in (:NamedTupleArray, :TupleArray, :AbstractTupleArray, :AbstractNamedTupleArray,
              :NamedTupleVector, :TupleVector)
    @eval _striptype(::Type{<:$_type}) = $_type
end

###
### DictViews interface
###

# DictViews.dictview_keys(ta::AbstractTupleArray) = first(tupgetdata(ta))

# function DictViews.dictview_values(ta::AbstractTupleArray)
#     n = length(tupgetdata(ta))
#     n == 1 && throw(ErrorException("there is no values field"))
#     n == 2 && return tupgetdata(ta)[2]
#     return tupgetdata(ta)[2:end]
# end

# # get a value by key
# # Default is slow linear search
# function DictViews.dictview_get(ta::AbstractTupleArray, k, default)
#     for i in eachindex(ta)
#         tup = ta[i]
#         if first(tup) == k
#             x = ta[i]
# #            return Base.rest(x, 2) # This is ok too, no faster
#             return ((x[i] for i in 2:length(x))...,)
#         end
#     end
#     return default
# end

end # module TupleArrays
