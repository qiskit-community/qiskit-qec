# Code circuit library
cmake_minimum_required(VERSION 3.12)
project(Circuits)

set(CIRCUITS_SRC
    intern/arctools.cpp
    
    arctools.h
)

add_library(libcircuits
    STATIC
    ${CIRCUITS_SRC}
)

target_include_directories(libcircuits
    PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

if (NOT (MSVC))
    target_compile_options(libcircuits PRIVATE -fno-strict-aliasing -fPIC ${ARCH_OPT})
else ()
    target_compile_options(libcircuits PRIVATE -fPIC ${ARCH_OPT})
endif ()

# Set prefix to "" since lib is already in project name
set_target_properties(libcircuits PROPERTIES PREFIX "")
