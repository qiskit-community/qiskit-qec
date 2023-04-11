"""
    module Bits

Tools for bits and masks. This module is similar to the registered Bits.jl

* Building masks (efficiently):
  `mask`, `rightmask`, `leftmask`, and `rangemask`.

* Bit views and indexing into bits of `Real`s:
  `bits`, `bit`

* Diagnostic: `min_bits`, `is_bitstring`, `bitsize`.

* Converting between representations of bit sequences:
  `bit_string`, `bool_tuple`, `bit_vector`
"""
module Bits

import Bits as _Bits
import Bits: bits, bitsize, BitVector1Mask, BitVector1, AbstractBitVector1
export bits, bitsize, BitVector1Mask, BitVector1, AbstractBitVector1
export undigits, undigits2

export mask, leftmask, rightmask, rangemask, asint, bit, min_bits

export bit_string, is_bitstring, bit_vector, bool_tuple

const Word = UInt128

"""
    rightmask([T=UInt128], i)

Return `n::T` such that the `i`th bit and all bits to the right (lower)
are one, and all bits to the left of the `i`th bit (higher) are zero.

See `leftmask`, `rangemask`, `mask`.
# Examples
```julia-repl
julia> bit_string(rightmask(UInt8, 3))
"00000111"
```
"""
rightmask(i) = rightmask(Word, i)
rightmask(::Type{T}, i::Integer) where {T} = one(T) << i - one(T)

"""
    leftmask([T=UInt128], i)

Return `n::T` such that the `i`th bit and all bits to the left (higher)
are one, and all bits to the right of the `i`th bit (lower) are zero.

See `rightmask`, `rangemask`, `mask`.
# Examples
```julia-repl
julia> bit_string(leftmask(UInt8, 3))
"11111100"
```
"""
leftmask(i) = leftmask(Word, i)
leftmask(::Type{T}, i::Integer) where {T} = ~(one(T) << (i - 1) - one(T))

"""
    rangemask([T=UInt128], ilo, ihi)
    rangemask([T=UInt128], (ilo, ihi)...)

Return `n::T` such that all bits in the range `ilo` to `ihi`, inclusive, are
one, and all other bits are zero. If `Tuples` `(ilo, ihi)...` are given,
then set bits in each range to one.

See `leftmask`, `rightmask`, `mask`.

# Examples
```julia-repl
julia> bit_string(rangemask(UInt8, 2, 7))
"01111110"

julia> bit_string(rangemask(UInt16, (1, 3), (5, 8), (14, 16)))
"1110000011110111"

julia> bit_string(rangemask(UInt8, (1, 5), (4, 8)))
"11111111"
```
"""
rangemask(args...) = rangemask(Word, args...)
rangemask(::Type{T}, ilo, ihi) where {T} = leftmask(T, ilo) & rightmask(T, ihi)
rangemask(::Type{T}, ranges::NTuple{2}...) where {T} =
    mapfoldl(x -> rangemask(T, x...), |, ranges)

"""
    mask([T=UInt128], i::Integer)
    mask([T=UInt128], r::UnitRange)
    mask([T=UInt128], itr)

Return `n::T` with the specified bits set to one and the rest zero.
The bits specified by each item in `itr` will be set to one.
Overlaps between ranges will have their bits set to one.

See `leftmask`, `rightmask`, `rangemask`.
# Examples
```julia-repl
julia> bit_string(mask(UInt8, 3))
"00000100"

julia> bit_string(mask(UInt8, (1, 5, 8)))
"10010001"

julia> bit_string(mask(UInt8, (2, 5, 8)))
"10010010"

julia> bit_string(mask(1:2:64))
"0101010101010101010101010101010101010101010101010101010101010101"

julia> bit_string(mask(UInt16, (1:3, 9, 14:16)))
"1110000100000111"
```
"""
mask(args) = mask(Word, args)
mask(::Type{T}, ur::UnitRange) where {T} = rangemask(T, ur.start, ur.stop)
# Following method has typical problem. For literal or::OR, it happens at compile time.
# For or::OR in variable it is 10x slower than the fallback below with `inds`.
# But, the fallback is slower for literal range, does not happen at compile time.
# mask(::Type{T}, or::OrdinalRange) where T = mask(T, collect(or))
mask(::Type{T}, i::Integer) where {T} = (one(T) << (i - one(T)))
mask(::Type{T}, inds) where {T} = mapfoldl(x -> mask(T, x), |, inds)
mask(::Type{T}, r::Base.OneTo) where {T} = rightmask(T, length(r))

# Modified from Bits.jl
asint(x::Integer) = x
asint(x::AbstractFloat) = reinterpret(Signed, x) # Signed gets all the types we need

"""
    bit(x::Real, i::Integer)

Similar to `Bits.bit` from registered `Bits.jl` package. A difference is that
the return type here is always `Int`.
"""
bit(x::Integer, i::Integer) = (x >>> UInt(i - 1)) & 1
bit(x::AbstractFloat, i::Integer) = bit(asint(x), i)
bit(x::Union{BigInt,BigFloat}, i::Integer) = Int(Base.GMP.MPZ.tstbit(x, i - 1))
bit(x::AbstractBitVector1, i::Integer) = x[i]
bit(x::AbstractVector{Bool}, i::Integer) = x[i]

"""
    bits(x::T, n) where {T<:Real}

Create a view of the bits of `x` that is truncated to the first `n` bits.

The underlying data is of type `T`. The truncated bits are zeroed.
Indexing at positions larger than `n` will not error, but will return zero.

# Examples
```julia-repl
julia> bits(0x5555)
<01010101 01010101>

julia> bits(0x5555, 10)
<01 01010101>

julia> bits(bits(0x5555, 10))
<00000001 01010101>
```
"""
_Bits.bits(x::T, n) where {T<:Real} = _Bits.BitVector1Mask(x & rightmask(T, n), n)

# TODO: Why do I do the following method? Why not reutrn x instead of x.x?
_Bits.bits(x::_Bits.BitVector1Mask) = _Bits.bits(x.x)
_Bits.bits(x::_Bits.BitVector1) = x

_Bits.bits(x::_Bits.BitVector1Mask, n) = _Bits.bits(x.x, n)
_Bits.bits(x::_Bits.BitVector1, n) = bits(x.x, n)

Base.count_ones(x::AbstractBitVector1) = Base.count_ones(x.x)

min_bits(x::_Bits.AbstractBitVector1) = min_bits(x.x)

# TODO: make this more robust. detect needed data type, etc.
"""
    bits([T=UInt128], s::AbstractString)

Convert the string `s` representing a binary number a one-word bit vector. This
method can be used to convert the string representation of
an `b::AbstractBitVector1` back to `b`. So, it behaves like `Base.parse(::Int,x)`
in that the bits in the string are in the same order as the bits in
`b`.

Spaces, and the characters '>', '<' are stripped from `s`.

# Examples
```julia-repl
julia> bits("00101")
<00101>

julia> bits("<01000000 11111111>")
<01000000 11111111>
```
"""
_Bits.bits(s::AbstractString) = _Bits.bits(UInt128, s)
function _Bits.bits(::Type{T}, s::AbstractString) where {T}
    ssimp = replace(s, r"[<>\s]" => "")
    x = parse(T, ssimp, base = 2)
    return _Bits.bits(x, length(ssimp))
end


const _VEC_LIKE = Union{
    AbstractVector{<:Integer},
    NTuple{<:Any,<:Integer},
    Base.Generator{<:AbstractVector},
}

"""
    bits([IntT=Int], dts::Union{AbstractVector{<:Integer},  NTuple{<:Any, <:Integer}}, n=length(dts))

Convert the container of binary digits `dts` to a `BitsVector1Mask`.

`IntT` is the storage type, i.e., the type that is wrapped. Input is not validated for
correctness, nor for having length greater than the number of bits in `IntT`.

# Examples
```julia-repl
julia> bits((0,1,0,1))
<0101>
```
"""
_Bits.bits(_digits::_VEC_LIKE, n = length(_digits)) = bits(Int, _digits, n)
_Bits.bits(::Type{IntT}, _digits::_VEC_LIKE, n = length(_digits)) where {IntT} =
    bits(undigits(IntT, _digits; base = 2), n)

"""
    undigits([IntT=Int], A; base=10)

The inverse of `digits`. That is `undigits(digits(n; base); base) == n`.
The number returned is of type `IntT` provided that the type is stable
under `+` and `*`, and that each element of `A` can be converted to `IntT`.
"""
undigits(A; base = 10) = undigits(Int, A, base = base)
function undigits(::Type{IntT}, A; base = 10) where {IntT}
    n = zero(IntT)
    @inbounds for i in reverse(eachindex(A))
        n = base * n + IntT(A[i])
    end
    return n
end

Base.zero(b::_Bits.AbstractBitVector1) = bits(zero(b.x))
Base.one(b::_Bits.AbstractBitVector1) = bits(one(b.x))

const _int_types =
    ((Symbol(pref, :Int, n) for n in (8, 16, 32, 64, 128) for pref in ("", "U"))...,)
for T in (_int_types..., :BigInt)
    @eval (Base.$T)(x::_Bits.AbstractBitVector1) = ($T)(x.x)
end
Base.Integer(x::_Bits.AbstractBitVector1) = x.x

for op in (:xor, :(&), :(|), :(+), :(-), :(*)) # it is actually useful sometimes to do +,-
    @eval (Base.$op)(x::_Bits.BitVector1Mask{T}, y::Real) where {T} =
        bits(($op)(x.x, T(y)), x.len)
    @eval (Base.$op)(x::_Bits.BitVector1{T}, y::Real) where {T} = bits(($op)(x.x, T(y)))
    @eval (Base.$op)(y::Real, x::_Bits.BitVector1{T}) where {T} = bits(($op)(x, y))
    @eval (Base.$op)(y::Real, x::_Bits.BitVector1Mask{T}) where {T} = bits(($op)(x, y))
    @eval (Base.$op)(y::BitVector1Mask, x::BitVector1Mask) = (
        (
            x.len == y.len ||
            throw(DimensionMismatch(string($op) * ": $(x.len) != $(y.len)"))
        );
        bits(($op)(x.x, y.x), x.len)
    )
end

# TODO: Is there and accessor method for .x ?
# We don't want lexicographic sorting. Also it uses slow fallback methods.
Base.:(<)(v1::_Bits.AbstractBitVector1, v2::_Bits.AbstractBitVector1) = v1.x < v2.x
Base.isless(v1::_Bits.AbstractBitVector1, v2::_Bits.AbstractBitVector1) = v1 < v2

# Could maybe use args... above for these
for op in (:(~),)
    @eval (Base.$op)(x::_Bits.BitVector1Mask{T}) where {T} = bits(($op)(x.x), x.len)
    @eval (Base.$op)(x::_Bits.BitVector1{T}) where {T} = bits(($op)(x.x))
end

Base.zero(b::_Bits.BitVector1Mask) = bits(zero(b.x), b.len)
Base.one(b::_Bits.BitVector1Mask) = bits(one(b.x), b.len)


function Base.reverse(b::BitVector1Mask{T}) where {T}
    c = zero(T)
    n = length(b)
    for i = 1:n
        c |= bit(b, i) << (n - i)
    end
    return bits(c, length(b))
end

###
### More bit functions
###

## Add more methods for Base.uinttype.
## Methods are defined in base for floating point types,
## and in the package BitIntegers.jl
let
    tups = [
        (Symbol(:Int, n), (Symbol(:UInt, n), Symbol(:Int, n), Symbol(:Float, n))) for
        n in (16, 32, 64)
    ]
    for (t, ts) in
        ((:Int8, (:UInt8, :Int8, :Bool)), tups..., (:Int128, (:UInt128, :Int128)))
        for tp in ts
            @eval inttype(::Type{$tp}) = $t
            @eval uinttype(::Type{$tp}) = $(Symbol("U", t))
        end
    end
end

"""
    bit_string(n::Integer; pad=nothing)
    bit_string(v; pad=0)

Give the literal bitstring representation of the number `n`. `bit_string` is
the same as `Base.bitstring`, except that the former allows specifying
left padding.

For `v` an iterable of known length, elements `x` in `v` are converted
characters `'0'` and `'1'` using `zero(x)` and `one(x)`.

`pad=0` ommits leading zeros in the output string

See also, `bit_vector`, `bool_vector`, and `bool_tuple`.

# Examples:

```julia-repl
julia> bit_string(128)
"0000000000000000000000000000000000000000000000000000000010000000"

julia> bit_string([1, 0, 1])
"101"

julia> bit_string(128; pad = 0)
"10000000"

julia> bit_string(128; pad = 9)
"010000000"
```
"""
function bit_string(x::T; pad::Union{Nothing,Integer} = nothing) where {T<:Integer}
    isnothing(pad) && return bitstring(x)
    return string(convert(uinttype(T), x); pad = pad, base = 2)
end

# TODO: make the semantics of `pad` the same as above for Integer input
function bit_string(v; pad::Union{Nothing,Integer} = 0)
    n = length(v)
    n_pad = pad - n > 0 ? pad - n : 0
    n_buf = n_pad > 0 ? pad : n
    buf = Vector{UInt8}(undef, n_buf)
    (zero_code, one_code) = (48, 49)
    @inbounds for i = 1:n_pad
        buf[i] = zero_code
    end
    j::Int = n_pad + 1
    @inbounds for i in eachindex(v)
        x = v[i]
        buf[j] =
            iszero(x) ? zero_code :
            isone(x) ? one_code : throw(DomainError(x, "Must be 0 or 1."))
        j += 1
    end
    return String(take!(IOBuffer(buf)))
end

"""
    min_bits(n::Integer)
    min_bits(bit_str::AbstractString)
    min_bits(v::Tuple)

Return the required number of bits in the binary representation
of `n` (or `bit_str`, or `v`).

The returned value is the position of the leftmost bit counting from the right equal to `1`.
Equivalently, the value is the length of the `bit_str` discounting leading zeros.

# Example
```julia-repl
julia> min_bits(2^10)
11

julia> min_bits("1"^10)
10

julia> min_bits([0,1,0])
2
```
"""
min_bits(n::T) where {T<:Integer} = 8 * sizeof(T) - Base.leading_zeros(n)

function min_bits(bit_str::AbstractString)
    n = Base.length(bit_str)
    for (i, c) in enumerate(bit_str)
        c === '1' && return (n - i + 1)
        i += 1
    end
    return 0
end

# could allow other iterables as well
function min_bits(v::Tuple)
    n = length(v)
    c::Int = 0
    for i in eachindex(v)
        @inbounds v[i] == one(v[i]) && break
        c += 1
    end
    return n - c
end

"""
    is_bitstring(bit_str::AbstractString; throw=false)

Return `true` if all characters in `bit_str` are either `'0'` or `'1'`,
otherwise `false`.
If `throw` is `true`, then throw an error rather than returning `false`.
"""
function is_bitstring(bit_str::AbstractString; throw = false)
    for c in bit_str
        if !(c === '1' || c === '0')
            throw && Base.throw(
                DomainError(c, "'$c' is not a bit. Characters must be one of ('0', '1')."),
            )
            return false
        end
    end
    return true
end

"""
    bool_tuple([IntT=Bool], bit_str::AbstractString)

Parse `bit_str` to a `Tuple` of `Bool`s.

If `IntT` is supplied then return a `Tuple` of `one(IntT)` and `zero(IntT)`.
`bit_str` is first validated.

# Examples

```julia-repl
julia> bool_tuple("10010")
(true, false, false, true, false)

julia> bool_tuple(Int, "10010")
(1, 0, 0, 1, 0)
```
"""
bool_tuple(bit_str::AbstractString) = bool_tuple(Bool, bit_str::AbstractString)
@inline function bool_tuple(::Type{IntT}, bit_str::AbstractString) where {IntT}
    is_bitstring(bit_str; throw = true)
    return ((c == '1' ? one(IntT) : zero(IntT) for c in bit_str)...,)
end

"""
    bit_vector(bit_str::AbstractString)

Parse `bit_str` to a `BitVector`. `bit_str` is first validated.
"""
bit_vector(bit_str::AbstractString) = BitVector(bool_tuple(bit_str))

"""
    bool_vector([IntT=Bool], bit_str::AbstractString)

Parse `bit_str` to a `Vector{Bool}`.

Return instead a `Vector{IntT}` if `IntT` is passed.
`bit_str` is first validated.
"""
bool_vector(bit_str::AbstractString) = bool_vector(Bool, bit_str)
function bool_vector(::Type{IntT}, bit_str::AbstractString) where {IntT}
    tup = bool_tuple(bit_str)
    vec = Vector{IntT}(undef, length(tup))
    for i = 1:length(tup) # eachindex ?
        @inbounds vec[i] = tup[i]
    end
    return vec
end

function bitvector1mask_2_string(bitvect1mask)
    stringbv1m = ""

    for val in bitvect1mask
        if val
            stringbv1m = "1" * stringbv1m
        else
            stringbv1m = "0" * stringbv1m
        end
    end
    return stringbv1m
end

end # module Bits
