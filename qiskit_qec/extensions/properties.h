#ifndef __Properties__
#define __Properties__

#include <vector>
#include <exception>

#include "linear.h"
#include "combinations.h"
#include "productiterator.h"

int minimum_distance(std::vector<std::vector<int> > &symplectic_vectors,
                     std::vector<std::vector<int> > &symplectic_gauge_vectors,
                     int max_weight=10);

#endif
