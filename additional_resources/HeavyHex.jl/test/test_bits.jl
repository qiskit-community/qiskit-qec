
using HeavyHex.Bits: min_bits, undigits, bits, bitvector1mask_2_string, BitVector1Mask

# This tests our extensions to Bits.jl

@testset "bits.jl" begin
    n = 12345
    @test undigits(digits(n)) == n
    @test undigits(digits(n, base = 2), base = 2) == n
    @test undigits(digits(n, base = 2), base = 10) != n
    tup = (1, 1, 1, 0, 0, 0)
    @test Tuple(bits(tup)) == tup
    @test Tuple(bits("000111")) == tup
    @test Tuple(bits("<000111>")) == tup
    @test string(bits(Tuple(bits("111000")))) == "<111000>"
    @test reverse(bits("111000")) == bits("000111")
    # 112 =1110000
    @test bitvector1mask_2_string(BitVector1Mask(112, 5)) == "10000"

end
