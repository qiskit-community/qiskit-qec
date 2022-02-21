#include "combinations.h"

Combinations::Combinations(int n_, int k_)
{
  if(n_ <= 0) throw std::invalid_argument("n must be positive");
  if(k_ < 0 || k_ > n_) throw std::invalid_argument("0 <= k <= n");
  n = n_;
  k = k_;
}

/*
 * Compute the next combination in lexicographic order.
 * Pass an empty vector "state" to begin.
 * Resume at any point by passing the current combination "state".
 * Behavior is undefined if "state" is not a valid lexicographically sorted combination.
 * Returns false when there are no more combinations.
 */
bool Combinations::next_combination(std::vector<int> &state)
{
  if(state.size() == 0)
  {
    for(int j=0; j<k; j++)
      state.push_back(j);
  } else {
    if(state[0] == n - k) return false;
    int i = k - 1;
    while(state[i] >= n - k + i) i--;
    state[i]++;
    for(int j=i+1; j<k; j++) state[j] = state[i] + j - i;
  }
  return true;
}

long int Combinations::size(void)
{
  int k0 = k;
  if(2*k > n) k0 = n - k;
  long int result = n;
  for(int i=2; i<=k0; i++)
  {
    result *= (n-i+1);
    result /= i;
  }
  return result;
}
