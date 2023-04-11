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
const MaybeNamedTuple{TupT} = Union{TupT,NamedTuple{<:Any,TupT}} where {TupT<:Tuple}

include("tuplearrays.jl")
include("bits.jl")
include("edges.jl")
include("io.jl")
include("experiment_helper.jl")
include("experiment.jl")
include("logging_utils.jl")
include("oml_decoding.jl")
include("job_runner.jl")
end
