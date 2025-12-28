# Copyright Matt Overby 2021.
# Distributed under the MIT License.

# Build visual/rendering examples
if (MCLGEOM_EXAMPLES)

    # Libigl for rendering/interface
    include(libigl)
    igl_include(glfw)

    add_executable(test_lbfgs ${CMAKE_CURRENT_SOURCE_DIR}/test_lbfgs.cpp)
    target_link_libraries(test_lbfgs PUBLIC igl::glfw mcl::geom)

endif()

# unit testing
if (MCLGEOM_UNIT_TESTS)

    enable_testing()

    find_package(TBB REQUIRED)

    add_executable(test_projection ${CMAKE_CURRENT_SOURCE_DIR}/test_projection.cpp)
    target_link_libraries(test_projection PUBLIC mcl::geom)
    add_test(NAME TestProjection COMMAND test_projection)

    add_executable(test_gaussseidel ${CMAKE_CURRENT_SOURCE_DIR}/test_gaussseidel.cpp)
    target_link_libraries(test_gaussseidel PUBLIC mcl::geom TBB::tbb)
    add_test(NAME TestGaussSeidel COMMAND test_gaussseidel)

    add_executable(test_graphcolor ${CMAKE_CURRENT_SOURCE_DIR}/test_graphcolor.cpp)
    target_link_libraries(test_graphcolor PUBLIC mcl::geom TBB::tbb)
    add_test(NAME TestGraphColor COMMAND test_graphcolor)

    add_executable(test_disjointsets ${CMAKE_CURRENT_SOURCE_DIR}/test_disjointsets.cpp)
    target_link_libraries(test_disjointsets PUBLIC mcl::geom)
    add_test(NAME TestDisjointSets COMMAND test_disjointsets)

endif()
