#include "arctools.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_c_circuits, module)
{
  module.doc() = "qiskit-qec code circuit extensions";
  module.def("_c_check_nodes", &check_nodes, "check_nodes in C++");
  module.def("_c_is_cluster_neutral", &is_cluster_neutral, "is_cluster_neutral in C++");
}