#include "faultenumerator.h"

  py::class_<FaultEnumerator>(m, "FaultEnumerator")
      .def(py::init<int, int, int, std::vector<std::vector<int>> &,
                    std::vector<int> &, std::vector<std::string> &,
                    std::vector<std::vector<std::string>> &>())
      .def("enumerate", &FaultEnumerator::enumerate)
      .def("reset", &FaultEnumerator::reset)
      .def("get_index", &FaultEnumerator::get_index)
      .def("get_state", &FaultEnumerator::get_state)
      .def("done", &FaultEnumerator::done);
