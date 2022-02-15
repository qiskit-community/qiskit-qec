#include "errorpropagator.h"
#include "faultsampler.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(compiledextension, m)
{
  m.doc() = "qiskit-qec compiled extensions";
  py::class_<ErrorPropagator>(m, "ErrorPropagator")

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
  py::class_<FaultSampler>(m, "FaultSampler")
      .def(py::init<int, int, std::vector<std::vector<int> >&,
                    std::vector<int>&, std::vector<std::string>&,
                    LabelToPauliWeightMap&,
                    std::map<std::string, double>&,
                    unsigned int>())
      .def("sample", &FaultSampler::sample);
}
