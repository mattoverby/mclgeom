# Copyright Matt Overby 2021.
# Distributed under the MIT License.
cmake_minimum_required(VERSION 3.0)

project(mclgeom_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

add_definitions(-DMCL_GEOM_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
option(MCL_BUILD_TESTS "Compile test" ON)

if (MCL_BUILD_TESTS)

	# Eigen
	find_package(EIGEN REQUIRED)
	include_directories(SYSTEM ${EIGEN_INCLUDE_DIR})

    # Should probably replace with actual unit testing through CTest
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
    add_executable(test_mcl ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.cpp)

endif()
