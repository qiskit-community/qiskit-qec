#ifndef __LinearSparse__
#define __LinearSparse__

#include <vector>
#include <exception>
#include <fstream> // Include for file I/O
#include <tuple>

/* actual interface */
std::vector<std::vector<int>> nullspace_sparse(const std::vector<std::vector<int>> &a, int n);
int nullity_sparse(const std::vector<std::vector<int>> &a, int n);
int rank_sparse(const std::vector<std::vector<int>> &a, int n);
bool lse_invalid_sparse(const std::vector<std::vector<int>> &a, const std::vector<bool> &b);
std::tuple<bool, std::vector<std::vector<int>>, std::vector<bool>> gaussian_elimination_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b);
std::vector<int> back_substitution_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b);
std::tuple<bool, std::vector<int>> solve_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b);
std::tuple<std::vector<int>, int, bool> minimize_weight_sparse(const std::vector<int> &x_part, const std::vector<std::vector<int>> &nullspace, int max_dof);
std::tuple<bool, std::vector<int>, int, bool> solve_optimal_sparse(std::vector<std::vector<int>> &a, int n, std::vector<bool> &b, int max_dof);
/* interface that is probably not needed, for debugging now */
std::vector<int> lor_sparse(const std::vector<int> &a1, const std::vector<int> &a2);
std::vector<int> land_sparse(const std::vector<int> &a1, const std::vector<int> &a2);
std::vector<int> lxor_sparse(const std::vector<int> &a1, const std::vector<int> &a2);
std::vector<int> lxor_sparse_a(std::vector<std::vector<int>> &a);
int lxor_weight_sparse(const std::vector<int> &a1, const std::vector<int> &a2);
int lxor_weight_sparse_a(std::vector<std::vector<int>> &a);
std::vector<int> ldif_sparse(const std::vector<int> &a1, const std::vector<int> &a2);
bool is_orthogonal_sparse(const std::vector<int> &a1, const std::vector<int> &a2);
bool in_support_sparse(const std::vector<int> &a, int x);

#endif
