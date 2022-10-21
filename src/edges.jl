module Edges

using ..TupleArrays: NamedTupleVector
using Dictionaries: Dictionary

export HEdges, error_set_prob

struct Probability{Td}
    prob_none::Float64 # prob that no error occurs
    corrections::Vector{Float64} # factors for each error that does occur
    dict::Td
end

function Probability(probs, bitstrings)
    prob_none = prod(1 - x for x in probs)
    corrections = [p / (1 - p) for p in probs]
#    dict = Dictionary(bitstrings, probs)
    dict = Dictionary(bitstrings, 1:length(bitstrings))
    return Probability(prob_none, corrections, dict)
end

struct HEdges{Tpr,Tb,Tpa,Td}
    probs::Tpr
    bitstrings::Tb
    pairs::Tpa
    pr_compute::Probability{Td}
end

function HEdges(eprobs::NamedTupleVector)
    bitstrings = eprobs(:bitstr)
    probs = eprobs(:prob)
    pr_compute = Probability(probs, bitstrings)
    return HEdges(probs, bitstrings, eprobs, pr_compute)
end

error_set_prob(edges::HEdges, error_strings) =
    error_set_prob(edges.pr_compute, error_strings)


function error_set_prob(pr::Probability, error_strings)
    p = pr.prob_none
    for estr = error_strings
#    for i in error_strings
        p1 = pr.dict[estr]
#        i = pr.dict[estr]
        p *= p1 / (1 - p1)
#        p *= pr.corrections[i]
    end
    return p
end


end #module Edges
