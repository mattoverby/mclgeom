// Copyright Matt Overby 2025.
// Distributed under the MIT License.

#ifndef MCL_GEOM_MULTICOLORGS_HPP
#define MCL_GEOM_MULTICOLORGS_HPP 1

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <functional>
#include <numeric>
#include <tbb/parallel_for.h>

namespace mcl {

template<typename T>
class MultiColorGaussSeidel
{
  public:
    using RowSparseMatrix = Eigen::SparseMatrix<T, Eigen::RowMajor>;

    /// @brief Solver options
    struct Options
    {
        int max_iters = 200; ///< max solver iters
        int check_tol = 10;  ///< num iters to check tol
        T rel_tol = (1e-8);  ///< delta x tolerance
        T omega = T(1.9);    ///< relaxation parameter
    };

    /// @brief Placeholder for post-sweep projection
    struct NoOp
    {
        template<typename DerivedX>
        void operator()(int, DerivedX&) const
        {
        }
    };

    /// @brief Solve AX=B s.t. X in C
    /// @param A sparse matrix
    /// @param B linear system rhs (row-major recommended)
    /// @param X linear system variable (row-major recommended)
    /// @param colors vertex colors for parallel execution
    /// @param parallel_exec optional flag for parallel/serial color
    /// @param project optional post-sweep projection (e.g., constraints)
    /// @param options solver settings
    /// @return number of iterations taken by solver
    template<typename DerivedX, typename DerivedB, typename Projection = NoOp>
    static int solve(const RowSparseMatrix& A,
                     const DerivedB& B,
                     DerivedX& X,
                     const std::vector<std::vector<int>>& colors = {},
                     const std::vector<bool>& parallel_exec = {},
                     const Projection& project = {},
                     const Options& options = Options())
    {
        int cols = X.cols();
        int num_colors = colors.size();
        int num_parallel_exec = parallel_exec.size();

        // Resize X to correct size if needed
        if (X.rows() != B.rows() || X.cols() != B.cols()) {
            X = Eigen::PlainObjectBase<DerivedX>::Zero(B.rows(), B.cols());
        }

        // If no colors, serial execution
        std::vector<int> all_inds;
        if (num_colors == 0) {
            num_colors = 1;
            all_inds.resize(X.rows());
            std::iota(all_inds.begin(), all_inds.end(), 0);
        }

        int iter = 1;
        for (; iter <= options.max_iters; ++iter) {
            for (int color = 0; color < num_colors; ++color) {

                const auto& inds = colors.empty() ? all_inds : colors[color];

                // Use serial execution if 1) no colors, 2) user-choice
                if (colors.empty() || (color < num_parallel_exec && !parallel_exec[color])) {
                    for (int ind : inds) {
                        sweep(ind, A, B, X, options);
                        project(ind, X);
                    }
                } else {
                    tbb::parallel_for(tbb::blocked_range<int>(0, int(inds.size())),
                                      [&](const tbb::blocked_range<int>& range) {
                                          for (int i = range.begin(); i != range.end(); ++i) {
                                              sweep(inds[i], A, B, X, options);
                                              project(inds[i], X);
                                          }
                                      });
                }
            }

            // Check if converged
            if (iter % options.check_tol == 0 && converged(A, B, X, options)) {
                break;
            }
        }
        return iter;
    }

    /// @brief Checks: ||b-Ax||/||b|| < rel_tol
    template<typename DerivedX, typename DerivedB>
    static bool converged(const RowSparseMatrix& A, const DerivedB& B, const DerivedX& X, const Options& options)
    {
        T B_norm = B.norm();
        T residual = (B - A * X).norm();
        return (residual / B_norm < options.rel_tol);
    }

    /// @brief Performs an X update at specified row
    template<typename DerivedX, typename DerivedB>
    static void sweep(int row, const RowSparseMatrix& A, const DerivedB& B, DerivedX& X, const Options& options)
    {
        T Aii = T(0);
        using RowVectorType = typename Eigen::internal::plain_row_type<DerivedX>::type;
        RowVectorType LUx = RowVectorType::Zero();
        for (typename RowSparseMatrix::InnerIterator iter(A, row); iter; ++iter) {
            if (iter.col() == row) {
                Aii = iter.value();
            } else {
                LUx += iter.value() * X.row(iter.col());
            }
        }
        if (std::abs(Aii) > T(0)) {
            auto X0 = X.row(row);
            auto X1 = (B.row(row) - LUx) / Aii;
            X.row(row) = (T(1) - options.omega) * X0 + options.omega * X1;
        }
    }
};

} // end namespace mcl

#endif
