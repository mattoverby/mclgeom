# Copyright Matt Overby 2021.
# Distributed under the MIT License.

cmake_minimum_required(VERSION 3.16)
project(mclgeom LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option(MCLGEOM_EXAMPLES "Compile visual examples" OFF)
option(MCLGEOM_UNIT_TESTS "Compile unit tests" OFF)

find_package(Eigen3 REQUIRED)

# Header only library
add_library(${PROJECT_NAME} INTERFACE)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
add_library(mcl::geom ALIAS ${PROJECT_NAME})

# testing
add_subdirectory(test)
