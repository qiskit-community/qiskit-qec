#include "errorpropagator.h"
#include "faultenumerator.h"
#include "faultsampler.h"
#include "distance.h"
#include "linear.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(_c_analysis, module)
{
  module.doc() = "qiskit-qec core analysis extensions";

  py::class_<ErrorPropagator>(module, "_CErrorPropagator")
      .def(py::init<int, int>())
      .def("load_circuit", &ErrorPropagator::load_circuit)
      .def("apply_error", &ErrorPropagator::apply_error)
      .def("cx", &ErrorPropagator::cx)
      .def("h", &ErrorPropagator::h)
      .def("s", &ErrorPropagator::s)
      .def("reset", &ErrorPropagator::reset)
      .def("measure", &ErrorPropagator::measure)
      .def("propagate", &ErrorPropagator::propagate)
      .def("get_cbits", &ErrorPropagator::get_cbits)
      .def("get_qubits", &ErrorPropagator::get_qubits)
      .def("get_qubit_array", &ErrorPropagator::get_qubit_array)
      .def("get_qreg_size", &ErrorPropagator::get_qreg_size)
      .def("get_creg_size", &ErrorPropagator::get_creg_size)
      .def("get_circuit_size", &ErrorPropagator::get_circuit_size)
      .def("__repr__",
           [](const ErrorPropagator &a)
           {
             return std::string("<ErrorPropagator with ") + std::to_string(a.get_qreg_size()) + std::string(" qubits, ") + std::to_string(a.get_creg_size()) + std::string(" bits, ") + std::to_string(a.get_circuit_size()) + std::string(" instructions>");
           });

  py::class_<FaultEnumerator>(module, "_CFaultEnumerator")
      .def(py::init<int, int, int, std::vector<std::vector<int>> &,
                    std::vector<int> &, std::vector<std::string> &,
                    std::vector<std::vector<std::string>> &>())
      .def("enumerate", &FaultEnumerator::enumerate)
      .def("reset", &FaultEnumerator::reset)
      .def("get_index", &FaultEnumerator::get_index)
      .def("get_state", &FaultEnumerator::get_state)
      .def("done", &FaultEnumerator::done);

  py::class_<FaultSampler>(module, "_CFaultSampler")
      .def(py::init<int, int, std::vector<std::vector<int>> &,
                    std::vector<int> &, std::vector<std::string> &,
                    LabelToPauliWeightMap &,
                    std::map<std::string, double> &,
                    unsigned int>())
      .def("sample", &FaultSampler::sample);

  module.def("_c_minimum_distance", &minimum_distance, "compute minimum distance of subsystem stabilizer code",
             py::arg("symplectic_vectors"), py::arg("symplectic_gauge_vectors"), py::arg("max_weight") = 10);
  module.def("_c_minimum_distance_by_tests", &minimum_distance_by_tests, "compute minimum distance of subsystem stabilizer code",
             py::arg("symplectic_vectors"), py::arg("symplectic_xl"), py::arg("symplectic_zl"), py::arg("max_weight") = 10);
  module.def("_c_rank", &rank, "compute rank of set of vectors", py::arg("vectors"));
  module.def("_c_isotropic", &is_isotropic, "test if set of symplectic vectors is isotropic", py::arg("symplectic_vectors"));
}
