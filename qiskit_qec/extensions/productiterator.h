#ifndef __ProductIterator__
#define __ProductIterator__

#include <vector>
#include <stdexcept>

// Iterate over a product of finite sets {0, 1, ..., n_i-1}
// for i = 1, 2, ..., m.
//
// Example: m = 2, n_1 = 4, n_2 = 3
// 00, 01, 02, 10, 11, 12, 20, 21, 22, 30, 31, 32
class ProductIterator
{
  private:
    std::vector<int> ni;

  public:
    ProductIterator(std::vector<int> &ni_);
    bool next_element(std::vector<int> &state);
    long int size(void);
};

#endif
