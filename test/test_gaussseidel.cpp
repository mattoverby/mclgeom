// Copyright Matt Overby 2021.
// Distributed under the MIT License.
#include <iostream>

#include <MCL/AssertHandler.hpp>
#include <MCL/GraphColor.hpp>
#include <MCL/MultiColorGaussSeidel.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <random>
#include <vector>

using RowSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using RowMatrixXd = Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;

int
gauss_seidel_no_colors();

int
gauss_seidel_colors();

int
main(int argc, char* argv[])
{
    (void)(argc);
    (void)(argv);
    int result_1 = gauss_seidel_no_colors();
    int result_2 = gauss_seidel_colors();
    return (result_1 == EXIT_SUCCESS && result_2 == EXIT_SUCCESS);
}

std::pair<RowSparseMatrix, RowMatrixXd> generate_linear_system(int n);

int
gauss_seidel_no_colors()
{
    int n = 10;
    auto [A, B] = generate_linear_system(n);

    // Solution
    Eigen::SimplicialLDLT<RowSparseMatrix> solver;
    solver.compute(A);
    mclAssert(solver.info() == Eigen::Success);
    RowMatrixXd X0 = solver.solve(B);
    mclAssert(X0.allFinite());

    // Test with Gauss Seidel solver.
    RowMatrixXd X;
    mcl::MultiColorGaussSeidel<RowMatrixXd> mcgs;
    mcgs.options.check_tol = 1; // check every iteration
    int iters = mcgs.solve(A, B, X);
    mclAssert(X.allFinite());
    mclAssert((X0 - X).lpNorm<Eigen::Infinity>() < 1e-6);

    // Starting closer to solution converges faster.
    RowMatrixXd X2 = X * 0.5;
    int iters2 = mcgs.solve(A, B, X);
    mclAssert(iters2 < iters);

    // Test with vector instead of matrix
    mcl::MultiColorGaussSeidel<Eigen::VectorXd> mcgs_vec;
    Eigen::VectorXd B_vec = B.col(0);
    Eigen::VectorXd X_vec;
    mcgs_vec.solve(A, B_vec, X_vec);
    mclAssert(X_vec.allFinite());
    mclAssert((X0.col(0) - X_vec).lpNorm<Eigen::Infinity>() < 1e-6);



#if 0
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
    mcl::MicroTimer t;
    graph_color.color();
    std::cout << "Color graph: " << t.elapsed_ms() << std::endl;

    std::vector<std::vector<int>> vertex_colors;
    graph_color.get_colors(vertex_colors);

    // Perform coloring with adjaceny matrix
    std::vector<std::vector<int>> vertex_colors_A;
    Eigen::SparseMatrix<double> D(T.rows(), V0.rows());
    D.setFromTriplets(adj_triplets.begin(), adj_triplets.end());
    Eigen::SparseMatrix<double> A = D.transpose() * D;
    t.reset();
    mcl::graph_color(A, vertex_colors_A);
    std::cout << "Make and color graph from matrix: " << t.elapsed_ms() << std::endl;

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

    check_colors(vertex_colors);
    check_colors(vertex_colors_A);

    // Combine small colors
    size_t min_color_size = 200;
    size_t prev_num_colors = vertex_colors.size();
    bool did_combine = mcl::combine_small_colors(min_color_size, vertex_colors);
    mclAssert(did_combine);
    mclAssert(vertex_colors.size() < prev_num_colors);

    // All colors smaller than min_color_size should be combined into
    // one color. That color might be less than the min size.
    for (int i = 0; i < int(vertex_colors.size() - 1); ++i) {
        mclAssert(vertex_colors[i].size() > min_color_size);
    }
#endif
    return EXIT_SUCCESS;
}

int
gauss_seidel_colors()
{
    return EXIT_SUCCESS;
}

std::pair<RowSparseMatrix, RowMatrixXd> generate_linear_system(int n)
{
    // 1D Poisson problem
    RowSparseMatrix A(n, n);
    std::vector<Eigen::Triplet<double>> triplets;
        for (int i = 0; i < n; ++i) {
            triplets.emplace_back(i, i, 4.0);
        if (i > 0) {
            triplets.emplace_back(i, i - 1, -1.0);
        }
        if (i < n - 1) {
            triplets.emplace_back(i, i + 1, -1.0);
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    RowMatrixXd B = RowMatrixXd::Ones(n, 3);
    return std::make_pair<RowSparseMatrix, RowMatrixXd>(std::move(A), std::move(B));
}