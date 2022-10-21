"""
    AnalyzeDistribution

These are tools that I thought at one time were part of the decoding
process. Turns out they are not. But, the results may still be of
some interest
"""
module AnalyzeDistribution

using Dictionaries: Dictionary
using Printf
using Combinatorics
import ..Bits
using ..Bits._Bits: BitVector1Mask, AbstractBitVector1, bits
using ..Edges: HEdges, error_set_prob

export checkerrors, checkcomb, write_example, read_example, analyze, analyze_all, print_report, analyze_analysis

"""
    write_example(fname, outcome_samp, errors_samp)

Write the example errors and observed outcome to a file.
"""
function write_example(fname, outcome_samp, errors_samp)
    open(fname, "w") do fh
        println(fh, "# observed outcome")
        println(fh, outcome_samp)
        println(fh)
        println(fh, "# secret errors")
        for s in errors_samp
            println(fh, s)
        end
    end
end

function read_example(fname)
    goodlines = []
    open(fname, "r") do fh
        while !eof(fh)
            line = readline(fh)
            match(r"\s*#", line) === nothing || continue
            match(r"^\s*$", line) === nothing || continue
            push!(goodlines, line)
        end
    end
    outcome_samp = bits(goodlines[1])
    errors_samp = bits.(goodlines[2:end])
    return (outcome_samp, errors_samp)
end

# prettier printing, like at the repl
replshow(x) = (show(stdout, MIME"text/plain"(), x); println())

function analyze_all(fname::AbstractString, edges::HEdges; num_errs_limit=4, info=true)
    (outcome_samp, errors_samp) = read_example(fname)
    return analyze_all(outcome_samp, errors_samp, edges; num_errs_limit, info)
end

function analyze_all(outcome_samp, errors_samp, edges; num_errs_limit=4, info=true)
    sample_prob = error_set_prob(edges, errors_samp)
    search_time = @elapsed prob_list = analyze(outcome_samp, edges; num_errs_limit, info)
    ind_found = search_correct_error_set(errors_samp, prob_list)
    num_errors_searched = length(filtererrors(outcome_samp, edges.pairs)) # TODO: optimize what is passed
    return (; outcome_samp, errors_samp, sample_prob, ind_found, prob_list, num_errors_searched,
            search_time, num_errs_limit) # named tuple
end

function analyze_analysis(analysis_result)
    (outcome_samp, errors_samp, sample_prob, ind_found, prob_list, num_errors_searched,
     search_time, num_errs_limit) = analysis_result
    none_found = isempty(prob_list)
    sets_too_small = (num_errs_limit < length(errors_samp))
    return (; ind_found, sets_too_small, none_found)
end


function print_report(analysis_result; print_limit=5)
    (outcome_samp, errors_samp, sample_prob, ind_found, prob_list, num_errors_searched,
     search_time, num_errs_limit) = analysis_result
    println("Observed bitstring:")
    println(outcome_samp, "\n")
    println("Number of error strings that have `1` at the same index (at least one) as the outcome string: $num_errors_searched\n")
    maxprob = isempty(prob_list) ? 1.0 : prob_list[1][:prob]
    println("Actual error string set, with probability = ", @sprintf("%.2e", sample_prob),
            ", relprob = ", @sprintf("%.2g", sample_prob/maxprob), " (relative to best candidate).")
    replshow(errors_samp)
    println()
    println("Error sets up to length $num_errs_limit searched in ", @sprintf("%.2g", search_time), " s\n")
    if length(errors_samp) > num_errs_limit
        printstyled("The actual error set could not be found because we did not search for large enough combinations.\n";
                    color=:read, bold=true)
    end
    goodcolor = :green # dark background
    if ind_found < 0
        printstyled("Actual error set not found in list of predictions.\n\n"; color=:red, bold=:true)
    else
        printstyled("Actual error *found* at prediciton number ", color=goodcolor)
        printstyled(ind_found; color=goodcolor, bold=true)
        printstyled(" (1 is best).\n\n"; color=goodcolor)
    end
    print_all_poss(prob_list; limit=print_limit)
end

"""
    analyze(outcome::AbstractBitVector1, edges::HEdges; num_errs_limit=4)

Attempt to find subsets of the set of error strings `edges` (and their probabilities) that
give `outcome` when `xor`d. That is when successively applied to flip bits. Search subsets
of errors up to size `num_errors_limit`. Four is usually a couple of seconds or less. Five
can take a couple of minutes.

The algorithm, which currently lacks finesse, is as follows.

* Find the subset, say `edge_pool` of `edges` that satisfy this criterion: There is at least one bit index `i` such that the error string has `1` at index `i` and that `outcome` also has a `1` at index `i`.

* Find all subsets of `edge_pool` of length `n` (with `n` taking successively values `1,2,3...`). `xor` each set and compare to `outcome`. That is, collect subsets of `edge_pool` (of size up to `num_errs_limit`) that could yield the observed `outcome` with non-zero probability.

* Compute the probabilty that each edge set occurs. Order these candidate sets by decreasing proability.

* Choose the first (most probable) on the list. This is useful because frequently the best candidate differs from the acutal error set and is also more probable than the actual error set.

When printing a report on the analysis in `print_report` we compute several other things, including the probability that the observed sample would occur.

# Notes

Of course, some error strings that are not in `edge_pool` may occur in the actual error
set. However, in a several tests, no such situation has produced a good candidate that is more
likely than the observed error set.  These error strings not in `edge_pool` are in two classes, some
flip may flip a bit that is `1` in the `outcome`.  Others only flip bits that are not `1` in
`outcome`. Cleary there is no algorithm that can tell whether they occured because `outcome` would
be the same if they were not present. Furthermore adding these errors to any set `X` does not change
`outcome`, but makes a set that is less probable than `X`.
"""
function analyze(outcome::AbstractBitVector1, edges::HEdges; num_errs_limit=4, info=true)
    prob_errors = filtererrors(outcome, edges.pairs) # TODO: optimize what is passed
#    errprobs_dicty = edges.dict
    combos_n = fill!(Vector{Any}(undef, 10), [])
    if info
        println("Size of error set,  number of combos of this size")
        println()
    end
    for nerrors in 1:num_errs_limit
        combos = checkcomb(prob_errors, outcome, nerrors)
        info && @show (nerrors, length(combos))
        combos_n[nerrors] = combos
    end
    prob_list = []
    for combos in combos_n
        for c in combos
            if !isempty(c)
                p = error_set_prob(edges, c)
                push!(prob_list, (prob=p, combo=c))
            end
        end
    end
#    sort!(prob_list; by=(x)->x[1], rev=true)
    sort!(prob_list; by=first, rev=true)
    isempty(prob_list) && return prob_list
    maxprob = prob_list[1][:prob]
    prob_list = [(prob=x[:prob], relprob=x[:prob]/maxprob, combo=x[:combo]) for x in prob_list]
    return prob_list
end

"""
    search_correct_error_set(errors_samp, prob_list)

Compare the actual sampled error set to our predictions, returning
the index of the prediction that matches, or `-1` if none matches.
The predictions are sorted in order of higher probability. So returned
index of `1` means the algorithm correctly predicted the error set.
"""
function search_correct_error_set(errors_samp, prob_list)
    n = -1
    for (i, error_res) in enumerate(prob_list)
        if error_res[:combo] == errors_samp
            n = i
            break
        end
    end
    return n
end

function print_one_poss(error_res)
    println("prob= ", @sprintf("%.2e", error_res[:prob]),
            ", relprob= ", @sprintf("%.2e", error_res[:relprob]),
            ", nerrs= ", length(error_res[:combo]))
    replshow(error_res[:combo])
end

function print_all_poss(prob_list; limit=10)
    n = length(prob_list)
    println("$n possible error combinations reported.")
    println("Printing $(min(length(prob_list),limit)) most probable out of $n")
    println()
    for (i, error_res) in enumerate(prob_list)
        i > limit && break
        print_one_poss(error_res)
        println()
    end
end

filtererrors(bstr, start) = filtererrors(((bstr, start),))

# Return *all* error strings that touch `1`s in `observed`.
function filtererrors(observed::AbstractBitVector1, errprobs)
    res = filter(x -> !iszero(x & observed), errprobs(:bitstr))
    [res...]
end

function filtererrors(tests::Union{AbstractVector, Tuple}, errprobs)
    fres = filter(matchat(tests), errprobs(:bitstr))
    [fres...]
end

"""
    checkerrors(estrings, samp)

Return `true if applying the sequence of errors `estrings` (to all-zero initial state)
gives the result `samp`.

That is, return `true` if the sequence of errors `estring` produces the observed
errors `samp`.
"""
function checkerrors(estrings::AbstractVector{<:AbstractBitVector1}, samp)
    observed = foldl(xor, estrings)
    return observed == samp
end

"""
    checkcomb(estrings, samp, n)

Return a `Vector` of all combinations of `n` error bitstrings
that produce observed string `samp`.

# Example
Generate all subsets of length `3` of `estrings` and return
only those which give `samp` when applied sequentially to
the all-zero initial state.
```julia-repl
julia> checkcomb(estrings, samp, 3)
```
"""
function checkcomb(estrings, samp, n)
    good = [] # eltype(estrings)[]
    for comb in combinations(estrings, n)
        if checkerrors(comb, samp)
            push!(good, comb)
        end
    end
    return good
end

# This is more heuristic and less exhaustive. It is one particular example
function analysis1()
    # Read saved output of: (errors_samp, outcome_samp) = one_sample(eprobs, Val(:all));
    (outcome, errors) = read_example("ex1.txt") # errors are secret
    @show outcome
    println()
    # Print indices of all 1s in outcome
    ones_pos = findall(isone, outcome)
    println("ones_pos = ", ones_pos)
    println()
    # Find all error strings that match the pattern below: "11" starting at position 9
    poss1 = filtererrors("11", 9)
    replshow(poss1)
    println()
    # Find all error strings that touch both of the two remaing bits
    poss2 = filtererrors((("1", 2), ("1", 24)))
    replshow(poss2)
    println()
    # We will try the union of these two sets
    poss3 = union(poss1, poss2)
    # It is evident that no single error is our outcome.
    # So, try all combinations of two errors from the set poss3
    # Return any combinations that give the output ovserved
    try2 = checkcomb(poss3, outcome, 2)
    @show length(try2)
    println()
    # Good! try1 holds a solution! Just one
    replshow(try2)
    println()
    # Looks correct, too. Check it with the secret solution
    @show errors == first(try2)
    # It is correct! ^^^
    # We should check the probabilties, which we have.
    # Are there combinations of three errors ? No.
    try3 = checkcomb(poss3, outcome, 3)
    @show length(try3)
    # How about combinations of 4 errors? Yes there are three of them
    try4 = checkcomb(poss3, outcome, 4)
    @show length(try4)
    # But, these are less probable than what we found. Again, we should
    # look at the probablilities.
    # We didn't do an exhastive search. Let's do this now.
    # The errors are not random looking. There are patterns we could exploit.
    # We did this heuristically above. But here, we use brute force.
    # We get all error strings that touch each 1 and then take their
    # union. We could do this all at once, as well.
    q2 = filtererrors("1", 2);
    q9 = filtererrors("1", 9);
    q10 = filtererrors("1", 10);
    q24 = filtererrors("1", 24);
    print(length.([q2, q9, q10, q24]))
    # [14, 36, 36, 10]
    qall = union(q2, q9, q10, q24);
    sum(length, [q2, q9, q10, q24])
    # 96
    length(qall)
    #86
    length(checkcomb(qall, outcome_samp, 1))
    # 0
    length(checkcomb(qall, outcome_samp, 2))
    # 1
    length(checkcomb(qall, outcome_samp, 3))
    # 12
    length(checkcomb(qall, outcome_samp, 4))
    # 117
    length(checkcomb(qall, outcome_samp, 5))
    # 1112
    nothing
end

end # module AnalyzeDistribution
