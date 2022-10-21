module DecoderIO

using ILog2: ILog2
using Dictionaries: Dictionary, sortkeys!
import Pickle
import Formatting

using ..DictTools: _AbstractDict, _insert!
using ..SmallDicts: SmallDictNamed, small_dict
using ..TupleArrays: NamedTupleVector, NamedTupleArray
using ..Bits: min_bits, bits
using ..Edges: HEdges
import ..HeavyHex
using BitIntegers: UInt256, UInt512, UInt1024

export read_hyper_edges, print_probs, write_ints_vector
export read_ints_vector, preprocessed_full_path, read_preprocessed, get_preprocessed_filename,
    preprocessed_tuples_to_bits

const ROOT_DATA_DIR = Ref(joinpath(dirname(dirname(pathof(HeavyHex))), "..", "HHCData"))
const PREPROCESSED_DIR = Ref(joinpath(ROOT_DATA_DIR[], "ps_counts_GMM_v4", "preprocessed_for_julia"))

function get_preprocessed_filename(job_id::Integer)
    fname = "PSall_and_last_init_both_fav2_false_3_3_3_counts_$(job_id)_preprocessed_plain_types_only.pkl"
    return joinpath(PREPROCESSED_DIR[], fname)
end

const TEST_DATA_DIR =
    joinpath(dirname(dirname(pathof(HeavyHex))),
        "examples", "test_data", "220107_Skywalker_data", "post_select_data",
        "ps_counts_GMM_v2")

preprocessed_full_path(fname) = joinpath(TEST_DATA_DIR, fname)

###
### DataFileSet
###

# Not flexible enough for templated directory, i.e. only filename is formatted.
struct DataFileSet
    root_dir::String
    set_dir::String
    file_template_fmt::Formatting.FormatExpr
end

function DataFileSet(root_dir::AbstractString, set_dir::AbstractString, file_template::AbstractString;
    rel_root=joinpath(dirname(dirname(pathof(HeavyHex))), ".."))
    if rel_root == nothing
        _root_dir = root_dir
    else
        _root_dir = joinpath(rel_root, root_dir)
    end
    return DataFileSet(_root_dir, set_dir, Formatting.FormatExpr(file_template))
end

# The default data file set
const DATA_FILE_SET = DataFileSet(
    "HHCData", joinpath("preprocessed", "plain_python_types"),
    "PSall_and_last_init_both_fav2_false_3_3_3_counts_{:d}_preprocessed_online.pkl"
)

# For example args can be a single integer argument: job_id
function get_preprocessed_filename(dfset::DataFileSet, args...)
    fname = Formatting.format(dfset.file_template_fmt, args...)
    return joinpath(dfset.root_dir, dfset.set_dir, fname)
end

get_preprocessed_filename(args...) = get_preprocessed_filename(DATA_FILE_SET, args...)

function read_preprocessed(dfset::DataFileSet, args...; format::Symbol=:auto, raw::Bool=false)
    fname = get_preprocessed_filename(dfset, args...)
    return read_preprocessed(fname; format=format, raw=raw)
end

read_preprocessed(job_id::Integer; format::Symbol=:auto, raw::Bool=false) =
    read_preprocessed(get_preprocessed_filename(job_id); format, raw)

function read_preprocessed(fname::AbstractString; format::Symbol=:auto, raw::Bool=false)
    if !(endswith(fname, ".pickle") || endswith(fname, ".pkl") || format === :pickle)
        throw(ErrorException("Unsupported format. If `fname` is a pickle file, use `format=:pickle`"))
    end
    data = Pickle.load(fname)
    data isa Dict || throw(ErrorException("Expecting a Dict from pickle file $fname"))

    (key1, key2) = ("preprocessed_info", "numerical_hyper_edges_dict")
    length(keys(data)) == 2 && collect(keys(data)) == [key1, key2] ||
        throw(ErrorException("Expecting two keys [\"$key1\", \"$key2\"]"))

    tup = (data[key1], data[key2])
    nt = NamedTuple{(Symbol(key1), Symbol(key2)),typeof(tup)}(tup)
    raw && return nt
    #    return preprocessed_tuples_to_bits(nt; ppinfo_type)
    return preprocessed_tuples_to_bits(nt)
end

function preprocessed_tuples_to_bits(pp::NamedTuple)
    num_bits = maximum(length(k) for k in keys(pp[1]))
    T = num_bits <= 128 ? UInt128 :
        error("Bit lengths greater than 128 are not supported.")

    return _preprocessed_tuples_to_bits(T, pp)
end

function _preprocessed_tuples_to_bits(::Type{T}, pp::NamedTuple) where {T}
    (ppinfo, edges) = (pp...,)
    delete!(ppinfo, "Total")
    ppinfo_bits = (bits(T, t) for t in keys(ppinfo))
    ppinfo_counts = values(ppinfo)

    new_ppinfo = NamedTupleVector((:bits, :counts), collect(ppinfo_bits), collect(ppinfo_counts))
    sort!(new_ppinfo)

    edges_bits = (bits(T, t) for t in keys(edges))
    probs = values(edges)
    new_edges = NamedTupleArray((:bits, :prob), collect(edges_bits), collect(probs))
    sort!(new_edges)

    new_tup = (new_ppinfo, new_edges)
    return NamedTuple{keys(pp),typeof(new_tup)}(new_tup)
end

# TODO: clean this up a bit Probably just return one of them and do conversion outside
# also separate conversion to bits
# function old_preprocessed_tuples_to_bits(pp::NamedTuple; ppinfo_type=:tuplevector, edges_type=:tuplevector)
#     (ppinfo, edges) = (pp...,)
#     delete!(ppinfo, "Total")
#     local new_ppinfo
#     if ppinfo_type === :tuplevector
#         new_ppinfo = NamedTupleVector((:bits, :counts), collect(ppinfo_bits), collect(ppinfo_counts))
#         sort!(new_ppinfo)
#     else
#         ppinfo_dictionary = sortkeys!(Dictionary(ppinfo_bits, ppinfo_counts))
#         if ppinfo_type === :vector
#             new_ppinfo = [(bits=t, count=n) for (t, n) in pairs(ppinfo_dictionary)]
#         elseif ppinfo_type in (:dict, :dictionary)
#             new_ppinfo = ppinfo_dictionary
#         else
#             error("unsupported type")
#         end
#     end
#     edges_bits = (bits(UInt, t) for (t, n) in pairs(edges))
#     probs = (p for (t, p) in pairs(edges))
#     if edges_type === :tuplevector
#         new_edges = NamedTupleArray((:bits, :prob), collect(edges_bits), collect(probs))
#         sort!(new_edges)
#     else
#         edge_dict = sortkeys!(Dictionary(edges_bits, probs))
#         if edges_type === :dict
#             new_edges = edge_dict
#         elseif edges_type === :vector
#             new_edges = [(bits=t, prob=p) for (t, p) in pairs(edges)]
#         else
#             error("unsupported type")
#         end
#     end
#     new_tup = (new_ppinfo, new_edges)
#     return NamedTuple{keys(pp),typeof(new_tup)}(new_tup)
# end


"""
    read_hyper_edges([IntT=Int64], fname)

Read hyperedge probabilities from a text file into a dictionary.

Items in the first column, the bit strings, are converted to `Int64` and are the
keys. The second column, the probabilities, are the values.  The columns are
whitespace separated.

Lines whose first non-whitespace character is `#` are ignored.

If the first argument is a type it will be used to represent bitstrings rather
than `Int64`.

If the first `1` in the bitstring is higher than the 63rd place (for `Int64`)
counting from the right, then the string cannot be represented
by an `Int64` and an error will be thrown.

The bitstrings were creaetd from Python tuples, with the bits ordered the same
in the string and the tuple. When converting a tuple to `AbstractBitsVector1`, our
fast bits representation, the bits are in reverse order. So, we reverse the string
after reading. This gives the same bit representation that we would get if we
read a tuple.
"""
read_hyper_edges(fname) = read_hyper_edges(UInt128, fname)
function read_hyper_edges(::Type{KeyT}, fname) where {KeyT}
    error_strings = Vector{KeyT}()
    probs = Vector{Float64}()
    open(fname, "r") do io
        for line in readlines(io)
            isnothing(match(r"\s*#", line)) || continue
            (bitstr, pstring) = split(line)
            # Reverse because this is a string. If it were a Tuple, no.
            bitstr = reverse(bitstr) # In python this is converted to a Tuple
            if KeyT <: AbstractString
                skey = bitstr
            else
                skey = parse(KeyT, bitstr; base=2) # assume KeyT is an integer
            end
            push!(error_strings, skey)
            push!(probs, parse(Float64, pstring))
            # _insert!(hyper_dict, skey, parse(Float64, pstring))
        end
    end
    return (error_strings, probs)
end

read_hyper_edges(::Type{SmallDictNamed}, fname) = read_hyper_edges(SmallDictNamed, UInt128, fname)
function read_hyper_edges(::Type{SmallDictNamed}, ::Type{KeyT}, fname) where {KeyT}
    (error_strings, probs) = read_hyper_edges(KeyT, fname)
    _min_bits = maximum(min_bits, error_strings)
    return small_dict(bits.(error_strings, _min_bits), probs; keyname=:bitstr, valuename=:prob)
end

read_hyper_edges(::Type{NamedTupleVector}, fname) = read_hyper_edges(NamedTupleVector, UInt128, fname)
function read_hyper_edges(::Type{NamedTupleVector}, ::Type{KeyT}, fname) where {KeyT}
    (error_strings, probs) = read_hyper_edges(KeyT, fname)
    _min_bits = maximum(min_bits, error_strings)
    return NamedTupleVector((:bitstr, :prob), (bits.(error_strings, _min_bits)...,), (probs...,))
end

function read_hyper_edges(::Type{HEdges}, fname::AbstractString)
    nt = read_hyper_edges(NamedTupleVector, fname)
    return HEdges(nt)
end

"""
    print_probs(prob_dist)

Print the iterable of probabilities `prob_dist` rendering the indices as bitstrings.
"""
print_probs(prob_dist::_AbstractDict, args...; pad=nothing) = _print_probs(prob_dist, args..., pad)
print_probs(prob_dist::AbstractVector, args...; pad=nothing) = _print_probs(pairs(prob_dist), args..., pad)

function print_probs(prob_dist::Vector{<:Pair}, args...; pad=nothing)
    _pad = isnothing(pad) ? ILog2.ilog2(maximum(first, prob_dist), RoundUp) : pad
    _print_probs(prob_dist, _pad)
end

function _print_probs(prob_dist, pad)
    _pad = isnothing(pad) ? ILog2.ilog2(maximum(keys(prob_dist)), RoundUp) : pad
    for (i, p) in prob_dist
        p > 0 && println(bit_string(i - 1; pad=_pad), ": ", p)
    end
end

"""
    write_ints_vector(fname::AbstractString, v::AbstractVector{Int})
    write_ints_vector(io::IO, v::AbstractVector{Int})

Write `v` in binary format to file `fname`.

Write the length of `v` as an `Int64`. Then write `v`. This is not
a preferred interface. But, JLD2 was updated during a session by
another Julia process and it broke JLD2 for that session. So, I
wrote this.
"""
function write_ints_vector(fname::AbstractString, v::AbstractVector{Int})
    open(fname, "w") do fh
        write_ints_vector(fh, v)
    end
end

function write_ints_vector(io::IO, v::AbstractVector{Int})
    write(io, length(v))
    write(io, v)
end

"""
    read_ints_vector(fname::AbstractString)

Read a vector of `Int64` from `fname`. One `Int`
will be read into `n`. Then a `Vector{Int}` of length
`n` will be read.
"""
function read_ints_vector(fname::AbstractString)
    local v
    open(fname, "r") do fh
        v = read_ints_vector(fh)
    end
    return v
end

function read_ints_vector(io::IO)
    n = read(io, Int)
    buf = Vector{Int}(undef, n)
    read!(io, buf)
    return buf
end

end  # module DecoderIO
