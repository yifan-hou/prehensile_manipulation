cmake_minimum_required(VERSION 2.8.12)
project(shared_grasping)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


add_compile_options(-std=c++11)

##
## Includes
##

find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)
find_library(CDDLIB cdd REQUIRED)



# qhull
SET(qhullDir /home/yifan/Git/qhull)
INCLUDE_DIRECTORIES(${qhullDir}/src/libqhullcpp)
INCLUDE_DIRECTORIES(${qhullDir}/src)
LINK_DIRECTORIES(${qhullDir}/build)
SET(qhullLibs qhullcpp qhull_r)





include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include/)


##
## targets
##
add_library(timer src/timer.cpp)
add_library(hfvc src/solvehfvc.cpp)

add_library(polylib src/polyhedron.cpp)
target_link_libraries(polylib ${CDDLIB} ${qhullLibs})

add_library(sharedgraspinglib src/wrench_space_analysis.cpp)
target_link_libraries(sharedgraspinglib polylib timer hfvc)

##
## Tests
##
add_executable(qhulltest src/testQHull.cpp)
target_link_libraries(qhulltest ${qhullLibs})

add_executable(polytest src/testpoly.cpp)
target_link_libraries(polytest polylib)

pybind11_add_module(example src/testPybind.cpp)
target_link_libraries(example PRIVATE sharedgraspinglib)

add_executable(lutest src/testEigenLU.cpp)
# add_executable(cddtest src/testcdd.cpp)
# target_link_libraries(cddtest ${CDDLIB})

