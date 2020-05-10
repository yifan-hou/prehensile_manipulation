cmake_minimum_required(VERSION 2.8.12)
project(shared_grasping)

add_compile_options(-std=c++11)

##
## Includes
##

find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)
find_library(CDDLIB cdd REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include/)


##
## targets
##
pybind11_add_module(example src/test.cpp)

add_library(polylib src/polyhedron.cpp)
target_link_libraries(polylib ${CDDLIB})



add_executable(cddtest src/testcdd.cpp)
target_link_libraries(cddtest ${CDDLIB})

add_executable(polytest src/testpoly.cpp)
target_link_libraries(polytest polylib)
