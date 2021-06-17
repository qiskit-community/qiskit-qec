# Synthesis of fault tolerant error correction circuits

**Problem**: *Given a code and a set of parameters, construct a fault-tolerant error correction circuit*

There are by now many ways to synthesize a fault-tolerant error correction circuit. One class of approaches measures each gauge or stabilizer generator using a separate block of ancillary qubits.

- bare circuits
- Shor circuits
- flag circuits
- ancilla-decoding circuits

Another class of approaches measures O(n-k) stabilizers in parallel.

- Steane EC (CSS codes)
- Knill EC

Recently there are also approaches that measure subsets of stabilizers in parallel.