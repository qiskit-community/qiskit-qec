#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <iostream>
#include <map>
#include <set>
#include <vector>

std::vector<int> check_nodes(
    std::vector<tuple<int, int, int, bool>> nodes, bool ignore_extra_logicals, bool minimal,
    map<tuple<int, int>, set<int>> cycle_dict,
    std::vector<tuple<int, int>> link_graph,
    map<int, vector<int>> link_neighbors,
    std::vector<int> z_logicals
    )