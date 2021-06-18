# Faults and error propagation

**Problem**: *Given a Clifford circuit, propagate all single faults to the end. (Scroll through the result in a GUI? Output as a list?) Same with all pairs of faults, or more?* (proposed by Ted Yoder)

**Detailed Description** Let $C=C_mC_{m-1}...C_1$ be a Clifford circuit on $n$ qubits. Let $Q$ be an $n$ qubit Pauli operator. Define $Q_0=Q$ and for $j=1,2,...,m$ let 
$$
Q_k=(C_kC_{k-1}...C_1)\cdot Q\cdot (C_kC_{k-1}...C_1)^\dagger
$$

**Problem**: *Given a stabilizer circuit, insert faults according to some model and propagate errors to the output qubits and measurements**

