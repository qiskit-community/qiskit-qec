#include "productiterator.h"

ProductIterator::ProductIterator(std::vector<int> &ni_)
{
  for(int i=0; i < ni_.size(); i++)
    if(ni_[i] <= 0) throw std::invalid_argument("each element must be positive");
  ni = ni_;
}

/*
 * Iterate over the product of sets.
 * Pass an empty vector "state" to begin.
 * Resume at any point by passing the current "state".
 * Behavior is undefined if "state" is not a valid element of the product.
 * Returns false when there are no more elements (element becomes invalid).
 */
bool ProductIterator::next_element(std::vector<int> &state)
{
  if(state.size() == 0)
  {
    state = std::vector<int>(ni.size(), 0);
  } else {
    for(int i=ni.size()-1; i>=0; i--)
    {
      state[i]++;
      if(state[i] != ni[i]) break;
      if(i!=0) state[i] = 0;
    }
    if(state[0] == ni[0]) return false;
  }
  return true;
}

long int ProductIterator::size(void)
{
  long int result = 1;
  for(auto i=ni.begin(); i!=ni.end(); i++) result *= *i;
  return result;
}
