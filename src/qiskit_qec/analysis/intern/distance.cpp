#include "distance.h"


/*
 * Compute the minimum distance of a stabilizer code.
 * The stabilizer corresponds to the space W spanned by a
 * set of symplectic vectors. The gauge group corresponds to
 * a space V spanned by another set of symplectic vectors
 * (possibly the same as W). The weight of a symplectic vector
 * is defined as the weight of the corresponding Pauli operator.
 * This algorithm enumerates all low weight vectors x in F_2^{2n}
 * to discover a lowest-weight vector in W^\perp - V.
 *
 * Returns the minimum distance of the corresponding code, or 0 if
 * the minimum distance exceeds max_weight.
 */
int minimum_distance(std::vector<std::vector<int> > &symplectic_vectors,
                     std::vector<std::vector<int> > &symplectic_gauge_vectors,
                     int max_weight) {
  if (!is_isotropic(symplectic_vectors))
    throw std::logic_error("symplectic vectors must span isotropic space");
  // We don't test: W is the center of V
  int n = symplectic_vectors[0].size() / 2;
  int rk = rank(symplectic_vectors);  // n - k - r
  int rk2 = rank(symplectic_gauge_vectors);  // n - k + r
  int k = n - (rk + rk2)/2;
  for (int weight=1; weight <= (max_weight < n ? max_weight : n); weight++) {
    std::vector<int> state;
    Combinations c(n, weight);
    std::vector<int> num_paulis(weight, 3);
    while (c.next_combination(state)) {
      std::vector<int> pstate;
      ProductIterator p(num_paulis);
      while (p.next_element(pstate)) {
        std::vector<int> error(2*n, 0);
        for (int i=0; i < weight; i++) {
          if (pstate[i] == 0 || pstate[i] == 2) error[state[i]] = 1;
          if (pstate[i] == 1 || pstate[i] == 2) error[state[i] + n] = 1;
        }
        // Test if error is in W^\perp
        bool inWperp = is_orthogonal(error, symplectic_vectors);
        // Test if error is in V
        std::vector<std::vector<int> > extended = symplectic_gauge_vectors;
        extended.push_back(error);
        bool inV = rank(extended) == rk2;
        if (inWperp)
          if ((k > 0 && !inV) || (k == 0 && inV))
            return weight;
      }
    }
  }
  return 0;
}

/*
 * Tests if a low weight logical Pauli operator exists.
 *
 * Tests whether there exists a Pauli operator with Hamming weight
 * <= weight that commutes with every row of the stabilizer table
 * (symplectic_vectors) and anti-commutes with a given logical operator
 * (symplectic_logical_op).
 *
 * Pauli operators on n qubits are represented as integer arrays
 * of length 2n in the order [X-part, Z-part].
 *
 * Returns True or False.
 */
bool distance_test(std::vector<std::vector<int> > &symplectic_vectors,
                   std::vector<int> &symplectic_logical_op,
                   int weight)
{
  int n = symplectic_vectors[0].size() / 2;  // num. qubits
  int m = symplectic_vectors.size();
  std::vector<int> pow2;
  for (int i = 0; i < m+1; i++)
    pow2.push_back(1 << i);
  int w1 = 0;
  int w2 = 0;
  if (weight % 2 == 0) {
    w1 = weight / 2;
    w2 = weight / 2;
  } else {
    w1 = (weight + 1) / 2;
    w2 = (weight - 1) / 2;
  }
  // Compute syndromes of all single-qubit errors
  std::vector<int> single_qubit_syndromes;
  for (int q=0; q < n; q++) {
    int syndX = 0;
    for (int i=0; i < m; i++) {
      if (symplectic_vectors[i][n + q] == 1)
        syndX += pow2[i];
    }
    if (symplectic_logical_op[n + q] == 1)
      syndX += pow2[m];
    int syndZ = 0;
    for (int i=0; i < m; i++) {
      if (symplectic_vectors[i][q] == 1)
        syndZ += pow2[i];
    }
    if (symplectic_logical_op[q] == 1)
      syndZ += pow2[m];
    single_qubit_syndromes.push_back(syndX);
    single_qubit_syndromes.push_back(syndZ);
    single_qubit_syndromes.push_back(syndX ^ syndZ);
  }
  // Examine all errors with weight w1
  int mask1 = 1 << m;
  std::set<int> T1c;
  std::set<int> T1a;
  if (w1 > 0) {
    Combinations c(single_qubit_syndromes.size(), w1);
    std::vector<int> state;
    while (c.next_combination(state)) {
      int synd = 0;
      for (int i=0; i < w1; i++)
        synd ^= single_qubit_syndromes[state[i]];
      if (synd & mask1)
        T1a.insert(synd ^ mask1);
      else
        T1c.insert(synd);
    }
  }
  // Examine all errors with weight w2
  std::set<int> T2c;
  std::set<int> T2a;
  if (w1 != w2) {
    if (w2 > 0) {
      Combinations c(single_qubit_syndromes.size(), w2);
      std::vector<int> state;
      while (c.next_combination(state)) {
        int synd = 0;
        for (int i=0; i < w2; i++)
          synd ^= single_qubit_syndromes[state[i]];
        if (synd & mask1)
          T2a.insert(synd ^ mask1);
        else
          T2c.insert(synd);
      }
    }
  } else {
    T2c = T1c;
    T2a = T1a;
  }
  std::set<int> intersect;
  std::set_intersection(T1c.begin(), T1c.end(), T2a.begin(), T2a.end(),
                        std::inserter(intersect, intersect.begin()));
  if (intersect.size() > 0) return true;
  intersect.clear();
  std::set_intersection(T1a.begin(), T1a.end(), T2c.begin(), T2c.end(),
                        std::inserter(intersect, intersect.begin()));
  if (intersect.size() > 0) return true;
  return false;
}

/*
 * Compute the minimum distance of a stabilizer code by
 * applying the distance test for each logical operator
 * in a basis.
 *
 * symplectic_vectors represent the stabilizer.
 * symplectic_xl represent the X logical operators.
 * symplectic_zl represent the Z logical operators.
 * 
 * Returns the minimum distance of the corresponding code, or 0 if
 * the minimum distance exceeds max_weight.
 */
int minimum_distance_by_tests(std::vector<std::vector<int> > &symplectic_vectors,
                              std::vector<std::vector<int> > &symplectic_xl,
                              std::vector<std::vector<int> > &symplectic_zl,
                              int max_weight) {
    int weight = max_weight + 1;
    for (int row=0; row < symplectic_xl.size(); row++) {
      for (int w=1; w < max_weight + 1; w++) {
        if (distance_test(symplectic_vectors, symplectic_xl[row], w)) {
          if (w < weight) weight = w;
          break;
        }
      }
    }
    for (int row=0; row < symplectic_zl.size(); row++) {
      for (int w=1; w < max_weight + 1; w++) {
        if (distance_test(symplectic_vectors, symplectic_zl[row], w)) {
          if (w < weight) weight = w;
          break;
        }
      }
    }
    if (weight < max_weight + 1)
      return weight;
    else
      return 0;
}
