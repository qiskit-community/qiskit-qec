# Faults and error propagation

**Problem**: *Given a Clifford circuit, propagate all single faults to the end. (Scroll through the result in a GUI? Output as a list?) Same with all pairs of faults, or more?* (proposed by Ted Yoder)

**Detailed Description**: Let $C=C_mC_{m-1}...C_1$ be a Clifford circuit on $n$ qubits. Let $Q$ be an $n$ qubit Pauli operator. Define $Q_0=Q$ and for $j=1,2,...,m$ let 
$$
Q_k=(C_kC_{k-1}...C_1)\cdot Q\cdot (C_kC_{k-1}...C_1)^\dagger.
$$
1. Given $C$ and $Q$ return the list $[Q_0,Q_1,...,Q_m]$.
2. Given $C$ and $Q$ show $Q_J$ visually on the curuit $C$ as you step through $j$. This could be done within the QEC GUI or with a Juypter notebook.

![image](https://user-images.githubusercontent.com/57962926/122604360-d5020e00-d043-11eb-8d8e-c4ef223f9031.png)

3. Provide methods to pass all singleton, pairs, etc. (or other defined combinations) through the circuit $C$ and gather totols, samples or statistics.

**Problem**: *Given a stabilizer circuit, insert faults according to some model and propagate errors to the output qubits and measurements**

