"""
    module DictTools

Some tools to support Dict and Dictionary (and their abstract supertypes.)

`_AbstractDict` is a `Union`. A few pirate methods are made so that their
APIs are closer.

Applications are `count_map`, `count_map!`, `update_map`, `update!`.
"""
module DictTools


export count_map, count_map!, update!, baretype, baretypeof, add_ncounts!, hist_to_dist

import Dictionaries
using Dictionaries: AbstractDictionary, Dictionary, gettoken, gettokenvalue, settokenvalue!

# Add a layer to handle both kinds of dictionaries

"""
    _AbstractDict{T, V}

Either an `AbstractDict` or an `AbstractDictionary`. A union type
"""
const _AbstractDict{T, V} = Union{AbstractDict{T,V}, AbstractDictionary{T,V}}

"""
    _Dict{T, V}

Either a `Dict` or a `Dictionary`. A union type
"""
const _Dict{T, V} = Union{Dict{T,V}, Dictionary{T,V}} # TODO: paramaterize

# Type piracy to give Dict the same construction interface as Dictionary
Base.Dict(inds, vals) = Dict(zip(inds, vals))

# TODO: Some examples would make this more clear
"""
    baretype(::Type{T})

Return the typename of `T`. This will fail for some input. In particular,
and in general, if `T` is a `UnionAll` type. However, more robust methods
are defined for `Dict` and `Dictionary`.
"""
baretype(::Type{T}) where {T} = isempty(T.parameters) ? T : T.name.wrapper
baretype(::Type{<:Dict}) = Dict
baretype(::Type{<:Dictionary}) = Dictionary

baretypeof(x) = baretype(typeof(x))

_insert!(d::AbstractDict, k, v) = (d[k] = v)
_insert!(d::Dictionaries.AbstractDictionary, k, v) = Dictionaries.insert!(d, k, v)

# These functions should be collected elsewhere.
# They may be actually collected elsewhere.

# The parameter `F` is actually necessary to force specialization
@inline function update!(dict::AbstractDictionary{T, V}, _key::T, func::F, default) where {F, T, V}
    (hasval, token) = gettoken(dict, _key)
    if hasval
        settokenvalue!(dict, token, func(gettokenvalue(dict, token)))
    else
        insert!(dict, _key, default)
    end
end

# Stolen from StatsBase. This is faster than naive method
@inline function update!(dict::Dict{T, V}, _key::T, func::F, default) where {F, T, V}
    index = Base.ht_keyindex2!(dict, _key)
    if index > 0
        @inbounds dict.vals[index] = func(dict.vals[index])
    else
        @inbounds Base._setindex!(dict, default, _key, -index)
    end
end

# This is a bit slower than both specialized methods
@inline function update!(dict::_AbstractDict{T,V}, _key::T, func::F, default) where {F, T, V}
    if haskey(dict, _key)
        dict[_key] = func(get(dict, _key))
    else
        dict[_key] = default
    end
end

function update!(dict::_AbstractDict, _keys, func::F, default) where F
    for _key in _keys
        update!(dict, _key, func, default)
    end
    return dict
end

update_map(_keys, func::F, default) where F = update_map(Dictionary, _keys, func, default)

update_map(::Type{DictT}, _keys, func::F, default) where {F, DictT} =
    update!(DictT{typeof(first(_keys)), typeof(default)}(), _keys, func, default)

"""
    count_map([::Type{T}=Dictionary], itr)

Return a dictionary of type `T` whose keys are elements of `itr`
and whose values count how many times each occurs.

`StatsBase.countmap` differs in that it has an optimization for
integers and that `itr` of indetermiate length is first collected
internally.
"""
count_map(args...) = update_map(args..., v -> (v + 1), 1)

"""
    count_map!(dict, itr)

Add items in `itr` to the count map `dict`. `itr` is any iterable,
including single iterable items, like numbers.
"""
@inline count_map!(dict::_AbstractDict{<:Any,V}, _key) where {V} =
    update!(dict, _key, v -> (v + one(V))::V, one(V))

function add_ncounts!(dict::_AbstractDict{<:Any,V}, _key, ncounts) where V
    update!(dict, _key, v -> (v + V(ncounts)), V(ncounts))
end

"""
    hist_to_dist(dict)

Convert the histogram `dict` to a probability distribution.

`dict` is a dictionary of counts. The output is a dictionary
with the same keys, and counts normalized to a probability distribution.
"""
hist_to_dist(dict) = baretypeof(dict)(keys(dict), values(dict) ./ sum(values(dict)))

end # module DictTools
