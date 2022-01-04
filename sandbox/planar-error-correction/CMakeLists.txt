cmake_minimum_required(VERSION 3.12...3.18)
project(compiledextension LANGUAGES CXX)

SET(SOURCE_DIR "src/plerco/pec_cpp_library")

include_directories(${SOURCE_DIR})

set(PYBIND11_FINDPYTHON ON) #requires CMAKE 3.12, alternatively run cmake -DPYBIND11_FINDPYTHON=ON ..
set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
add_subdirectory(extern/pybind11)

pybind11_add_module(compiledextension "${SOURCE_DIR}/CompiledExtensions.cpp" "${SOURCE_DIR}/ErrorPropagator.cpp" "${SOURCE_DIR}/FaultEnumerator.cpp" "${SOURCE_DIR}/Combinations.cpp" "${SOURCE_DIR}/ProductIterate.cpp")
