// Copyright Matt Overby 2021.
// Distributed under the MIT License.
#include <MCL/MultiColorGaussSeidel.hpp>

#include <MCL/AssertHandler.hpp>
#include <MCL/GraphColor.hpp>
#include <MCL/MicroTimer.hpp>
#include <MCL/ShapeFactory.hpp>

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <random>
#include <vector>

using RowSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using RowMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

void
build_graph_laplacian(int nx, int ny, RowSparseMatrix& A, RowMatrixXd& B, std::vector<std::vector<int>>& colors);

void
test_gauss_seidel(const RowSparseMatrix& A, const RowMatrixXd& B, const std::vector<std::vector<int>>& colors);

int
main(int argc, char* argv[])
{
    (void)(argc);
    (void)(argv);

    int n = 10;
    RowSparseMatrix A;
    RowMatrixXd B;
    std::vector<std::vector<int>> colors; // red-black example
    build_graph_laplacian(n, n, A, B, colors);
    printf("A = %d x %d\n", int(A.rows()), int(A.cols()));
    mclAssert(colors.size() > 1);
    mclAssert(mcl::verify_graph_colors(Eigen::SparseMatrix<double>(A), colors));

    mcl::MicroTimer t;
    test_gauss_seidel(A, B, {}); // no colors (serial)
    printf("Serial GS: %f ms\n", t.elapsed_ms());

    t.reset();
    test_gauss_seidel(A, B, colors); // colors (parallel)
    printf("Parallel GS: %f ms\n", t.elapsed_ms());

    return EXIT_SUCCESS;
}

void
test_gauss_seidel(const RowSparseMatrix& A, const RowMatrixXd& B, const std::vector<std::vector<int>>& colors)
{
    // Solution
    Eigen::SimplicialLDLT<RowSparseMatrix> solver;
    solver.compute(A);
    mclAssert(solver.info() == Eigen::Success);
    RowMatrixXd X0 = solver.solve(B);
    mclAssert(X0.allFinite());

    // Test with Gauss-Seidel solver.
    RowMatrixXd X;
    mcl::MultiColorGaussSeidel<RowMatrixXd> mcgs;
    double tol = mcgs.options.rel_tol * 10;
    mcgs.options.max_iters = 100;
    mcgs.options.check_tol = 1; // check every iteration
    int iters = mcgs.solve(A, B, X, colors);
    mclAssert(X.allFinite());
    mclAssert(iters < 50);
    mclAssert((X0 - X).lpNorm<Eigen::Infinity>() < tol);

    // Starting closer to solution converges faster.
    RowMatrixXd X2 = X * 0.5;
    int iters2 = mcgs.solve(A, B, X, colors);
    mclAssert(iters2 < iters);

    // Test with vector instead of matrix, most common implementations
    mcl::MultiColorGaussSeidel<Eigen::VectorXd> mcgs_vec;
    Eigen::VectorXd B_vec = B.col(0);
    Eigen::VectorXd X_vec;
    int iters3 = mcgs_vec.solve(A, B_vec, X_vec, colors);
    mclAssert(X_vec.allFinite());
    mclAssert((X0.col(0) - X_vec).lpNorm<Eigen::Infinity>() < tol);
    mclAssert(iters3 < 50);

    // Test with projection
    X.setZero();
    mcgs.project = [](int i, RowMatrixXd& X_) { X_(i, 2) = std::max(X_(i, 2), 1.00001); };
    int iters4 = mcgs.solve(A, B, X, colors);
    mclAssert((X.col(2).array() > 1.0).all());
}

void
build_graph_laplacian(int nx, int ny, RowSparseMatrix& A, RowMatrixXd& B, std::vector<std::vector<int>>& colors)
{
    const int n = nx * ny;

    A.resize(n, n);
    B.resize(n, 3);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(5 * n);

    auto id = [nx](int x, int y) { return x + y * nx; };

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const int i = id(x, y);
            int deg = 0;

            auto add_edge = [&](int nx, int ny) {
                triplets.emplace_back(i, id(nx, ny), -1.0);
                ++deg;
            };

            if (x > 0) {
                add_edge(x - 1, y);
            }
            if (x < nx - 1) {
                add_edge(x + 1, y);
            }
            if (y > 0) {
                add_edge(x, y - 1);
            }
            if (y < ny - 1) {
                add_edge(x, y + 1);
            }

            triplets.emplace_back(i, i, static_cast<double>(deg));
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    for (int i = 0; i < A.rows(); ++i) {
        A.coeffRef(i, i) += 1.0;
    }

    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const int i = id(x, y);
            B(i, 0) = 1.0;
            B(i, 1) = static_cast<double>(x);
            B(i, 2) = static_cast<double>(y);
        }
    }

    colors.resize(2);
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            int color = (x + y) % 2;
            colors[color].emplace_back(id(x, y));
        }
    }
}
