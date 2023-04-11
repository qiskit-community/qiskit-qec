"""
module DecoderIO

Handles the I/O needed to read in the preprocessed experiment data files stored 
in a Julia compatible format in the the HHCData repo. 
"""
module DecoderIO

using ILog2: ILog2
using Dictionaries: Dictionary, sortkeys!
import Pickle
import Formatting
using Printf
using ..TupleArrays: NamedTupleVector, NamedTupleArray
using ..Bits: min_bits, bits
using ..Edges: HEdges
import ..HeavyHex
using BitIntegers: UInt256, UInt512, UInt1024
using Dictionaries
export read_hyperedges, write_ints_vector
export read_ints_vector,
    preprocessed_full_path, read_preprocessed, preprocessed_tuples_to_bits

"""Assumes all the data is being pulled from a copy of the HHCData repo"""
const FILENAME_TEMPLATE = joinpath(
    "ps_counts_GMM_v4",
    "preprocessed_for_julia",
    "PSall_and_last_init_both_fav2_false_3_3_3_counts_%d_preprocessed_plain_types_only.pkl",
)

""" Global  dynamic variable storing the path to the HHCData"""
const HHC_DIRECTORY_PATH = Ref("")

function setup_hhc_directory_path(hhcdata_directory_path::AbstractString)
    global HHC_DIRECTORY_PATH[] = hhcdata_directory_path
end

function load_job_run_from_id(job_id)

    preprocessed_data = read_preprocessed(job_id)

    _hyperedges = preprocessed_data.numerical_hyperedges_dict
    hyperedges = Dictionaries.Dictionary(_hyperedges.bits, _hyperedges.prob)

    prep_info = preprocessed_data.preprocessed_info
    observed_counts = prep_info.counts

    observations_all = reverse.(prep_info.bits)

    function bit_from_start(observations_all)
        observations = [x[2:end] for x in observations_all]
        expected_bits = [x[1:1] for x in observations_all]
        return (observations, expected_bits)
    end

    (observations, expected_bits) = bit_from_start(observations_all)
    return hyperedges, observed_counts, observations, expected_bits

end

"""
    read_preprocessed(job_id::Integer)::String

Get the preprocessed data file for `job_id` from the default data directory. This value is set by the user at runtime.
"""

function read_preprocessed(job_id::Int64; format::Symbol=:auto, raw::Bool=false)
    fname = Printf.format(Printf.Format(FILENAME_TEMPLATE), job_id)
    file_path = joinpath(HHC_DIRECTORY_PATH[], fname)
    return read_preprocessed(file_path; format, raw)
end

function read_preprocessed(fname::AbstractString; format::Symbol=:auto, raw::Bool=false)
    if !(endswith(fname, ".pickle") || endswith(fname, ".pkl") || format === :pickle)
        throw(
            ErrorException(
                "Unsupported format. If `fname` is a pickle file, use `format=:pickle`",
            ),
        )
    end
    data = Pickle.load(fname)
    data isa Dict || throw(ErrorException("Expecting a Dict from pickle file $fname"))

    (key1, key2) = ("preprocessed_info", "numerical_hyper_edges_dict")
    length(keys(data)) == 2 && collect(keys(data)) == [key1, key2] || throw(
        ErrorException(
            "Expecting two keys [\"$key1\", \"$key2\"]. Found $(length(keys(data))) keys.",
        ),
    )

    tup = (data[key1], data[key2])
    # Correct the spelling of the second element of the NamedTuple
    nt = NamedTuple{(Symbol(key1), :numerical_hyperedges_dict),typeof(tup)}(tup)
    raw && return nt
    #    return preprocessed_tuples_to_bits(nt; ppinfo_type)
    return preprocessed_tuples_to_bits(nt)
end
function preprocessed_tuples_to_bits(pp::NamedTuple)
    num_bits = maximum(length(k) for k in keys(pp[1]))
    T = num_bits <= 128 ? UInt128 : error("Bit lengths greater than 128 are not supported.")

    return _preprocessed_tuples_to_bits(T, pp)
end

function _preprocessed_tuples_to_bits(::Type{T}, pp::NamedTuple) where {T}
    (ppinfo, edges) = (pp...,)
    delete!(ppinfo, "Total")
    ppinfo_bits = (bits(T, t) for t in keys(ppinfo))
    ppinfo_counts = values(ppinfo)

    new_ppinfo =
        NamedTupleVector((:bits, :counts), collect(ppinfo_bits), collect(ppinfo_counts))
    sort!(new_ppinfo)

    edges_bits = (bits(T, t) for t in keys(edges))
    probs = values(edges)
    new_edges = NamedTupleArray((:bits, :prob), collect(edges_bits), collect(probs))
    sort!(new_edges)

    new_tup = (new_ppinfo, new_edges)
    return NamedTuple{keys(pp),typeof(new_tup)}(new_tup)
end

"""
    read_hyperedges([IntT=Int64], fname)

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

The bitstrings were created from Python tuples, with the bits ordered the same
in the string and the tuple. When converting a tuple to `AbstractBitsVector1`, our
fast bits representation, the bits are in reverse order. So, we reverse the string
after reading. This gives the same bit representation that we would get if we
read a tuple.
"""
read_hyperedges(fname) = read_hyperedges(UInt128, fname)
function read_hyperedges(::Type{KeyT}, fname) where {KeyT}
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
        end
    end
    return (error_strings, probs)
end

read_hyperedges(::Type{NamedTupleVector}, fname) =
    read_hyperedges(NamedTupleVector, UInt128, fname)
function read_hyperedges(::Type{NamedTupleVector}, ::Type{KeyT}, fname) where {KeyT}
    (error_strings, probs) = read_hyperedges(KeyT, fname)
    _min_bits = maximum(min_bits, error_strings)
    return NamedTupleVector(
        (:bitstr, :prob),
        (bits.(error_strings, _min_bits)...,),
        (probs...,),
    )
end

function read_hyperedges(::Type{HEdges}, fname::AbstractString)
    nt = read_hyperedges(NamedTupleVector, fname)
    return HEdges(nt)
end

end  # module DecoderIO
