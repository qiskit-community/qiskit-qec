# Analysis library
cmake_minimum_required(VERSION 3.12)
project(Analysis)

set(ANALYSIS_SRC
    intern/combinations.cpp
    intern/distance.cpp
    intern/errorpropagator.cpp
    intern/faultenumerator.cpp
    intern/faultsampler.cpp
    intern/linear.cpp
    intern/productiterator.cpp
    
    combinations.h
    distance.h
    errorpropagator.h
    faultenumerator.h
    faultsampler.h
    linear.h
    productiterator.h
)

add_library(libanalysis
    STATIC
    ${ANALYSIS_SRC}
)

target_include_directories(libanalysis
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

if (NOT (MSVC))
    target_compile_options(libanalysis PRIVATE -fno-strict-aliasing -fPIC ${ARCH_OPT})
else ()
    target_compile_options(libanalysis PRIVATE -fPIC ${ARCH_OPT})
endif ()

# Set prefix to "" since lib is already in project name
set_target_properties(libanalysis PROPERTIES PREFIX "")
