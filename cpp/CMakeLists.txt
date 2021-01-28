cmake_minimum_required(VERSION 2.8.12)
project(shared_grasping)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# debug & profiling
set(CMAKE_BUILD_TYPE Debug)
# SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
# SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
# SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

add_compile_options(-std=c++14)

##
## Includes
##

find_package(Eigen3 REQUIRED)
find_package(pybind11 REQUIRED)
find_library(PPLLIB ppl REQUIRED)
find_library(GMPLIB gmp REQUIRED)
find_library(GLPKLIB glpk REQUIRED)
find_package(modus REQUIRED)
find_package(ortools CONFIG REQUIRED)

# qhull
SET(qhullDir /home/yifan/Git/qhull)
INCLUDE_DIRECTORIES(${qhullDir}/src/libqhullcpp)
INCLUDE_DIRECTORIES(${qhullDir}/src)
LINK_DIRECTORIES(${qhullDir}/build)
SET(qhullLibs qhullcpp qhull_r)

# Robot utilities
find_library(RobotUtilLib Utilities HINTS /usr/local/lib/RobotUtilities)

# yaml
include_directories(yaml-cpp)
find_library(YAML_LIBRARY yaml-cpp HINTS /usr/local/lib)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include/)


##
## targets
##
add_library(timer SHARED src/timer.cpp)
add_library(hfvc SHARED src/solvehfvc.cpp)

add_library(polylib SHARED src/polyhedron.cpp)
target_link_libraries(polylib ${PPLLIB} ${GMPLIB} ${qhullLibs} ${GLPKLIB} ${RobotUtilLib})

add_library(jacobianlib SHARED src/shared_grasping_jacobian.cpp)

add_library(sharedgraspinglib SHARED src/wrench_space_analysis.cpp)
target_link_libraries(sharedgraspinglib polylib timer hfvc ${RobotUtilLib}
	${YAML_LIBRARY})

##
## Tests
##
# add_executable(qhulltest src/testQHull.cpp)
# target_link_libraries(qhulltest ${qhullLibs})

add_executable(polytest src/testpoly.cpp)
target_link_libraries(polytest polylib)

pybind11_add_module(wrenchStampingLib src/wrenchStampingPybind.cpp)
target_link_libraries(wrenchStampingLib PRIVATE sharedgraspinglib
	polylib timer hfvc ${RobotUtilLib}
	${YAML_LIBRARY})

add_executable(lutest src/testEigenLU.cpp)

add_executable(ppltest src/testppl.cpp)
target_link_libraries(ppltest ${PPLLIB} ${GMPLIB})

add_executable(lptest src/lp.cpp)
target_link_libraries(lptest ${GLPKLIB})

add_executable(sharedGraspingCpptest src/test_wrenchstamping.cpp)
target_link_libraries(sharedGraspingCpptest
	jacobianlib
	sharedgraspinglib
	modus
	ortools::ortools)


# # for debug
# find_library(libcnpy cnpy HINTS /usr/local/lib)
# add_executable(exam src/exam.cpp)
# target_link_libraries(exam sharedgraspinglib ${libcnpy})