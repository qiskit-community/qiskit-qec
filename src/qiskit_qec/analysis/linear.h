#ifndef __Linear__
#define __Linear__

#include <vector>
#include <exception>
#include <fstream> // Include for file I/O
#include <tuple>

int symplectic_inner_product(std::vector<int> symplectic_a,
                             std::vector<int> symplectic_b);
bool is_isotropic(std::vector<std::vector<int> > &symplectic_vectors);
bool is_orthogonal(std::vector<int> x,
                   std::vector<std::vector<int> > &symplectic_vectors);
int rank(std::vector<std::vector<int> > vectors);

std::tuple<bool, std::vector<std::vector<bool>>, std::vector<bool>, std::vector<bool>> solve(std::vector<std::vector<bool>> &a, std::vector<bool> &b);
bool contradiction(const std::vector<std::vector<bool>> &a, const std::vector<bool> &b);
//std::tuple<bool, std::vector<std::vector<int>>, std::vector<bool>, std::vector<int>> solve_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b);
// bool contradiction_sparse(const std::vector<std::vector<int>> &a, const std::vector<bool> &b);
// bool contains(const std::vector<int> &a, int x);
// std::vector<int> sparse_xor(const std::vector<int> &a1, const std::vector<int> &a2);
// std::vector<int> sparse_dif(const std::vector<int> &a1, const std::vector<int> &a2);
// std::vector<bool> back_propagation_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b);

// std::vector<std::vector<int>> kernel(std::vector<std::vector<int>> a, int n);

#endif
