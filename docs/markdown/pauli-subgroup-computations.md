## Pauli Subgroup Computations

This set of requirements covers the various group computations for subgroups of the $n$-qubit Pauli Groups $\mathcal{P}_n$.

---

**Problem**: Given a subgroup $A$ of $\mathcal{P}_n$ find a minimal generating set for $A$. 

Note: As $\mathcal{P}_n$ is a $p$-group it follows by Burnside's basis theorem that all minimal
generating sets (by including) have the samecardinality, hence there is not ambiguity in the question.

---

**Problem**: Given a subgroup $A$ of $\mathcal{P}_n$ find a hyperbolic generating set for $A$. 
An optional argument should return a normalized version.

A set $Q=\\{A_1,A_2,...,A_k,A_{k+1},B_{k+1},...,A_{k+t},B_{k+t}\\}\in\mathcal{P}_n$ is a hyperbolic generating
set if 
- $Q$ is a minimal generating set for $\langle Q\rangle$
- $Z(Q)=\langle A_1,A_2,...,A_k \rangle$
- $A_i$ commutes with all generators except for $B_i$.
- $B_i$ commutes with all generators except for $A_i$.

Thus the standard generating set $\\{iI,X_1,Z_1,X_2,Z_2,...,X_n,Z_n\\}$ is a hyperbolic generating set for
$\mathcal{P}_n$.

A hyperbolic generating set $Q=\\{A_1,A_2,...,A_k,A_{k+1},B_{k+1},...,A_{k+t},B_{k+t}\\}\in\mathcal{P}_n$ 
is said to be normal if all generators have order $2$ except possibly for one element of the center.

A hyperbolic generating set can be found using a modified Gram-Schmidt process.  See for example 
-T. Brun, Min-Hsiu Hsieh arxiv.orf:1610.04013

---

**Problem**: Let $A$ be a subgroup of $\mathcal{P}_n$ with hyperbolic generating set  $Q=\\{A_1,A_2,...,A_k,A_{k+1},B_{k+1},...,A_{k+t},B_{k+t}\\}$. Extend $Q$ to a hyperbolic generating set
for $\mathcal{P}_n$ 

---

**Problem**: Given a subset $X$ of $\mathcal{P}_n$ find its Center, Normalizer, Centralizer, relative to $\mathcal{P}_n$.

---

**Problem**: Given a Gauge group find a minimal generating set for its bare logical operators

---

**Problem**: Let $A$ be a subgroup of $\mathcal{P}_n$. 1) Find all maximal subgroups of $A$. 2) Find all maximal gauge groups of $A$, find all maximal stabilizer subgroups of $A$. Construct generating sets (hyperbolic, or normal hyperbolic is required) for all maximal subgroups found.



