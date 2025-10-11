# Copyright Matt Overby 2021.
# Distributed under the MIT License.

cmake_minimum_required(VERSION 3.16)
project(mclgeom LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(MCLGEOM_TESTING_ONLY "Compile tests, no rendering" OFF)

find_package(Eigen3 REQUIRED)

# Header only library
add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
add_library(mcl::geom ALIAS ${PROJECT_NAME})

# Build visual/rendering examples
if (NOT MCLGEOM_TESTING_ONLY)

    # Libigl for rendering/interface
    include(libigl)
	igl_include(glfw)

	# Set define for loading objs
	add_definitions(-DMCLGEOM_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")

    add_executable(lbfgs ${CMAKE_CURRENT_SOURCE_DIR}/test/lbfgs.cpp)
    target_link_libraries(lbfgs PUBLIC igl::glfw mcl::geom)
endif()

# unit testing (work in progress)
enable_testing()
add_executable(test_projection ${CMAKE_CURRENT_SOURCE_DIR}/test/test_projection.cpp)
target_link_libraries(test_projection PUBLIC mcl::geom)
add_test(NAME TestProjection COMMAND test_projection)