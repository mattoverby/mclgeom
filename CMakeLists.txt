# Copyright Matt Overby 2021.
# Distributed under the MIT License.
cmake_minimum_required(VERSION 3.16)
project(mclgeom_project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
add_definitions(-DMCL_GEOM_ROOT_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
option(MCL_GEOM_BUILD_TESTS "Compile tests" OFF)

if (MCL_GEOM_BUILD_TESTS)

    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
    option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
    option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
    include(libigl)

    # Should probably replace with actual unit testing through CTest
    add_executable(test_mcl ${CMAKE_CURRENT_SOURCE_DIR}/tests/test.cpp)
    target_link_libraries(test_mcl PUBLIC Eigen3::Eigen)
    add_executable(lbfgs ${CMAKE_CURRENT_SOURCE_DIR}/tests/lbfgs.cpp)
    target_link_libraries(lbfgs PUBLIC Eigen3::Eigen igl::core igl::opengl_glfw)

endif()
