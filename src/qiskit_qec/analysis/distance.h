#ifndef __Distance__
#define __Distance__

#include <algorithm>
#include <vector>
#include <set>
#include <exception>
#include <iterator>

#include "linear.h"
#include "combinations.h"
#include "productiterator.h"

int minimum_distance(std::vector<std::vector<int> > &symplectic_vectors,
                     std::vector<std::vector<int> > &symplectic_gauge_vectors,
                     int max_weight=10);

bool distance_test(std::vector<std::vector<int> > &symplectic_vectors,
                   std::vector<int> &symplectic_logical_op,
                   int weight);

int minimum_distance_by_tests(std::vector<std::vector<int> > &symplectic_vectors,
                              std::vector<std::vector<int> > &symplectic_xl,
                              std::vector<std::vector<int> > &symplectic_zl,
                              int max_weight);

#endif
