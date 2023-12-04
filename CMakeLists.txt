cmake_minimum_required(VERSION 3.4...3.18)
project(msproteomics)

add_subdirectory(pybind11)
pybind11_add_module(python_example src/siteloc.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(python_example
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
