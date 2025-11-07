// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include <Eigen/Sparse>
#include <MCL/AssertHandler.hpp>
#include <MCL/GraphColor.hpp>
#include <MCL/ReadEleNode.hpp>

#include <iostream>
#include <unordered_map>
#include <vector>

int
graph_color_mesh();

int
main(int argc, char* argv[])
{
    (void)(argc);
    (void)(argv);
    return graph_color_mesh();
}

int
graph_color_mesh()
{
    using namespace Eigen;
    typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;
    typedef Matrix<int, Dynamic, Dynamic, RowMajor> RowMatrixXi;
    RowMatrixXd V0;
    RowMatrixXi T;

    // Load mesh
    {
        MatrixXd inV;
        MatrixXi inT;
        if (!mcl::read_ele_node(MCLGEOM_ROOT_DIR "/test/armadillo_3k", inV, inT)) {
            std::cout << "Failed to load " << MCLGEOM_ROOT_DIR "/test/armadillo_3k" << std::endl;
            return EXIT_FAILURE;
        }
        V0 = inV;
        T = inT;
    }

    mclAssert(T.rows() > 0 && T.cols() == 4);
    mclAssert(V0.rows() > 0 && V0.cols() == 3);
    mcl::GraphColor graph_color(V0.rows());

    // Create an adjacency matrix. Okay to have redundant elements.
    std::vector<Eigen::Triplet<double>> adj_triplets;
    for (int i = 0; i < T.rows(); ++i) {
        for (int j = 0; j < 4; ++j) {
            int t0 = T(i, j);
            int t1 = T(i, (j + 1) % 4);
            adj_triplets.emplace_back(i, t0, 1);
            adj_triplets.emplace_back(i, t1, 1);
            graph_color.make_union(t0, t1);
        }
    }

    // Perform coloring
    std::vector<std::vector<int>> vertex_colors;
    graph_color.color();
    graph_color.get_colors(vertex_colors);

    // Perform coloring with adjaceny matrix
    std::vector<std::vector<int>> vertex_colors_A;
    Eigen::SparseMatrix<double> D(T.rows(), V0.rows());
    D.setFromTriplets(adj_triplets.begin(), adj_triplets.end());
    Eigen::SparseMatrix<double> A = D.transpose() * D;
    mcl::graph_color(A, vertex_colors_A);

    auto check_colors = [&](const std::vector<std::vector<int>>& colors) {
        // For easier verification, move to map.
        Eigen::VectorXi vertex_assigned = Eigen::VectorXi::Zero(V0.rows());
        std::unordered_map<int, int> vertex_to_color;
        for (int i = 0; i < (int)colors.size(); ++i) {
            mclAssert(colors.size() > 0);
            for (int vtx : colors[i]) {
                mclAssert(vertex_to_color.count(vtx) == 0);
                vertex_to_color[vtx] = i;
                vertex_assigned[vtx] = 1;
            }
        }
        mclAssert((vertex_assigned.array() > 0).all());

        // Verify: no vertex in same color that is a part of the same tet stencil(s)
        for (int i = 0; i < T.rows(); ++i) {
            for (int j = 0; j < 4; ++j) {
                int c0 = vertex_to_color[T(i, j)];
                int c1 = vertex_to_color[T(i, (j + 1) % 4)];
                mclAssert(c0 != c1);
            }
        }
    };

    std::cout << "colors: " << vertex_colors.size() << ", " << vertex_colors_A.size() << std::endl;
    check_colors(vertex_colors);
    check_colors(vertex_colors_A);

    return EXIT_SUCCESS;
}