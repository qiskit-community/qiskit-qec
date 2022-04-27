#ifndef __Linear__
#define __Linear__

#include <vector>
#include <exception>

int symplectic_inner_product(std::vector<int> symplectic_a,
                             std::vector<int> symplectic_b);
bool is_isotropic(std::vector<std::vector<int> > &symplectic_vectors);
bool is_orthogonal(std::vector<int> x,
                   std::vector<std::vector<int> > &symplectic_vectors);
int rank(std::vector<std::vector<int> > vectors);

#endif
