"""
    module OmlDecoding

    Contains the Maximum Likelihood Decoder 
"""
module OmlDecoding

using Dictionaries: Dictionary, Dictionaries, set!, unset!, insert!, delete!
using ..Edges: HEdges
import ..DecoderIO
import ..Bits
import ..Experiment
using ConcurrentCollections
# export stuff just for dev to get completion at the repl
export OnlineMLDecoder, setup_job, update, decode, test_job_info

struct OnlineMLDecoder{T,F,K,V}
    numerical_hyperedges::T
    bit_order_func::F
    zero_distribution::K
    zero_error_key::Any
    prefix_length::Int
    logical_correction::Bool
    prefix_store::V
end

function OnlineMLDecoder(jobinfo::NamedTuple; hyperedges = nothing)
    local _hyperedges
    _hyperedges = hyperedges

    prefix_length = jobinfo.prefix_len
    logical_correction = true
    prefix_store = ConcurrentCollections.ConcurrentDict{keytype(_hyperedges),Any}()
    hyperedge_keys = collect(keys(_hyperedges))
    L = length(hyperedge_keys[1])
    k1 = Bits.bits(UInt128(0), L)
    zero_distribution = Dictionary{typeof(k1),Float64}([k1], [1.0])
    zero_error_key = k1


    return OnlineMLDecoder(
        hyperedges,
        jobinfo.bit_order_func,
        zero_distribution,
        zero_error_key,
        prefix_length,
        logical_correction,
        prefix_store,
    )
end

function setup_job(dec::OnlineMLDecoder)
    dec.prefix_store = ConcurrentCollections.ConcurrentDict{keytype(_hyperedges),Any}() # make new dict
end

function decode(decoder::OnlineMLDecoder, observation::Tbs) where {Tbs}
    len_obs = length(observation)
    # figure out the most advanced start in the prefix_dictionary
    _bit = 0
    prefix::Tbs = observation[[len_obs - decoder.bit_order_func(i) for i = 0:len_obs-1]]
    local _keys::Vector{Tbs}# assert type to catch errors
    local distribution::Dictionary{Tbs,Float64}
    for i = decoder.prefix_length:-1:1
        pref1i::Tbs = prefix[1:i]
        value = ConcurrentCollections.maybeget(decoder.prefix_store, pref1i)
        if value !== nothing
            (_keys, distribution) = something(value)
            _bit = i
            break
        end
    end
    if iszero(_bit)
        _keys = collect(keys(decoder.numerical_hyperedges))
        L = length(_keys[1]::Tbs)
        k1 = Bits.bits(UInt128(0), L)::Tbs
        distribution = Dictionary{typeof(k1),Float64}([k1], [1.0]) # {tuple(np.zeros(L)): 1}
    end
    jkeys = empty(_keys)
    for j = _bit:(len_obs-1)
        newkeys = empty(_keys) # needs to be new, not just a buffer. It is stored below
        # jkeys = empty(_keys)
        empty!(jkeys) # Not much performance gain here.
        thisbit = decoder.bit_order_func(j)
        for k in reverse(eachindex(_keys)) # (length(_keys)-1):-1:0
            if _keys[k][thisbit+1] == 1
                push!(jkeys, _keys[k])
            else
                push!(newkeys, _keys[k])
            end
        end
        _keys = newkeys
        # include self.numerical_hyperedges_dict
        # Some kind of aliasing going on between these dictionaries.
        # It's baffling.

        ############ surgery zone begin
        jdistribution =
            Dictionary(Dict(key => decoder.numerical_hyperedges[key] for key in jkeys)) # TODO this might be causing a big slowdown

        new_hyperedges_1, new_hyperedges_2 = split_hyperedges_dictionary(jdistribution)
        new_distribution =
            build_distribution_w_large_merge(new_hyperedges_1, new_hyperedges_2, decoder)
        distribution = process_merge_two_dicts(new_distribution, distribution) #will this be more accurate but slower ?

        ########### surgery zone end
        # expand the self.prefix_dictionary
        # truncate
        key_set = Set(keys(distribution))
        for k in key_set
            if k[thisbit+1] != prefix[j+1]
                delete!(distribution, k)
            end
        end

        # TODO if I don't write to prefix_store before I'm done messing  with distribution, this will probably stop trying to retrieve non-existent keys.
        if j < decoder.prefix_length
            # the value for any key is predetermined so it's fine (I think) if two separate threads both write to the same key as long as they do so atomically to avoid undefined behavior
            decoder.prefix_store[prefix[1:j+1]] = (_keys, distribution)
        end
    end
    if !decoder.logical_correction
        return distribution
    end
    if isempty(distribution)
        println("distribution is empty. Something is wrong.")
        return nothing
    end
    maxkey = findmax(distribution)[2] # findmax returns (v, k), with max v.
    return (maxkey[len_obs+1:end], distribution)
end

update(pdict, hyperedge) = update(pdict, typeof(pdict)(), hyperedge)
function update(pdict, newdict, hyperedge)
    (v, p) = hyperedge
    empty!(newdict)
    for k in keys(pdict)
        val = pdict[k] * (1 - p)
        insert!(newdict, k, val)
    end
    for k in keys(pdict)
        newk = xor(k, v) # tuple([(k[i] + v[i]) % 2 for i in range(len(k))])
        if newk in keys(newdict)
            Dictionaries.set!(newdict, newk, newdict[newk] + pdict[k] * p) # same thing ?
        else
            val = pdict[k] * p
            Dictionaries.set!(newdict, newk, val)
        end
    end
    return newdict
end


function split_hyperedges_dictionary(hyperedges)
    hyperedges_1 = Dictionary()
    hyperedges_2 = Dictionary()

    i = 0
    for (h, p) in pairs(hyperedges)
        if i == 1
            insert!(hyperedges_1, h, p)
            i = 0
        else
            insert!(hyperedges_2, h, p)
            i = 1
        end

    end
    return hyperedges_1, hyperedges_2
end

function build_distribution_via_linear_merge_from_hyperedges(
    hyperedges,
    decoder::OnlineMLDecoder,
)
    cur_hyperedge_distributions = [
        Dictionary(Dict(decoder.zero_error_key => 1 - p, h => p)) for
        (h, p) in pairs(hyperedges)
    ]
    og_distr = Dictionary(Dict(decoder.zero_error_key => 1.0))
    for hd in cur_hyperedge_distributions
        og_distr = process_merge_two_dicts(og_distr, hd)
    end
    return og_distr
end

function build_distribution_w_large_merge(h1, h2, decoder::OnlineMLDecoder)
    d1 = build_distribution_via_linear_merge_from_hyperedges(h1, decoder::OnlineMLDecoder)
    d2 = build_distribution_via_linear_merge_from_hyperedges(h2, decoder::OnlineMLDecoder)
    dtot = process_merge_two_dicts(d1, d2)
    return dtot
end

function process_merge_two_dicts(dict1::Dictionary, dict2::Dictionary)
    merged = typeof(dict1)()
    for (k1, v1) in pairs(dict1)
        for (k2, v2) in pairs(dict2)
            new_key = xor(k1, k2) #
            if !(new_key in keys(merged))
                insert!(merged, new_key, v1 * v2)
            else
                set!(merged, new_key, merged[new_key] += v1 * v2)
            end
        end
    end
    return merged
end



end # module OmlDecoding
