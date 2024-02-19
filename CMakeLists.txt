cmake_minimum_required(VERSION 3.12)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR
        "Building of the code within the source directory is not currently allowed"
        "\n Please remove any CMakeCache.txt files created in the root of the source directory"
    )
endif()

# ------------------------------------------------------------------------------
# Policies
#
# Make any need policy settings here. This should be updated once a final
# structure is approved i.e. use policies to inforce the rules!

cmake_policy(VERSION 3.12)

# ------------------------------------------------------------------------------
# Project Initialization

project(Qiskit_qec)

# ------------------------------------------------------------------------------
# Macros
#
# Load up any macros
# Project level macros should be place in build_files/cmake/macros.cmake

include(build_files/cmake/macros.cmake)

# ------------------------------------------------------------------------------
# Set Options

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------------------------
# Set compile flags

if (NOT(MSVC))
    set(ARCH_OPT "-O3" "-mno-avx2")
else ()
    set(ARCH_OPT "-O2")
    # https://stackoverflow.com/a/8591946 (RTC1 incompatible with O2 flag)
    STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif ()


# ------------------------------------------------------------------------------
# Add code subdirectories to make libraries

# Source CMakeLists
add_subdirectory(intern)
add_subdirectory(extern)
add_subdirectory(src)

# Make Pybind python bindings

set(ANALYSIS_BINDINGS_SRC "src/qiskit_qec/analysis/bindings/analysis_bindings.cpp")
pybind11_add_module(_c_analysis ${ANALYSIS_BINDINGS_SRC})
target_link_libraries(_c_analysis PRIVATE libanalysis)

set(CIRCUITS_BINDINGS_SRC "src/qiskit_qec/circuits/bindings/circuits_bindings.cpp")
pybind11_add_module(_c_circuits ${CIRCUITS_BINDINGS_SRC})
target_link_libraries(_c_circuits PRIVATE libcircuits)



