module HeavyHex

export MaybeNamedTuple

"""
    MaybeNamedTuple{TupT}

A `DataType` that may be either a `NamedTuple` or a `Tuple`.

# Examples
```julia-repl
julia> tups = [(1,2), (a=1, b=2), (1.0, 2), (a=1.0, b=2)];

julia> print(isa.(tups, MaybeNamedTuple))
Bool[1, 1, 1, 1]

julia> print(isa.(tups, MaybeNamedTuple{Tuple{Int, Int}}))
Bool[1, 1, 0, 0]
```
"""
const MaybeNamedTuple{TupT} = Union{TupT,NamedTuple{<:Any,TupT}} where TupT <: Tuple

# const IntegerLike = Union{Integer, AbstractBitVector1}

include("tuplearrays.jl")
include("dict_tools.jl")
include("small_dict.jl")
include("bits.jl")
include("utils.jl")
include("edges.jl")
include("io.jl")
include("analyzer_helper.jl")
include("preprocessor_helper.jl")
include("experiment.jl")
include("decode.jl")
include("analyze_distribution.jl")
include("bigdecoder.jl")
include("montecarlo.jl")
include("analyze.jl")

end
