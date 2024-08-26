#ifndef __Combinations__
#define __Combinations__

#include <vector>
#include <stdexcept>

// Enumerate subsets of {0, 1, 2, ..., n} of size k
// where each subset is sorted in increasing order
// and the subsets are ordered in lexicographic order.
//
// In lexicographical order, a subset a_1 a_2 ... a_k < b_1 b_2 ... b_k
// if a_i < b_i on the first index i where a_i and b_i differ.
//
// Example: n = 4, k = 2
// Choose 2 elements from {0, 1, 2, 3}.
// 01 < 02 < 03 < 12 < 13 < 23
class Combinations
{
  private:
    int n;
    int k;

  public:
    Combinations(int n_, int k_);
    bool next_combination(std::vector<int> &state);
    long int size(void);
};

#endif
