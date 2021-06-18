# Lookup table decoding

**Problem**: *Given a generating set of a stabilizer group, construct a lookup table decoder for Pauli noise*

For each Pauli error E of weight less than O(d), we record the pair (s(E), E') where s(E) is the syndrome of E and E' is the lowest weight error with that syndrome. We ignore errors F such that wt(F) > wt(E) and s(E)=s(F).

Given the automorphism group of a code, it may not be necessary to generate and store the entire lookup table. A special case occurs for cyclic codes called the Meggitt decoder.

This decoding method is not scalable or optimal but can be applied to any small code.