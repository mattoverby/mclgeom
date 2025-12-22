// Copyright Matt Overby 2025.
// Distributed under the MIT License.

#ifndef MCL_GEOM_MULTICOLORGS_HPP
#define MCL_GEOM_MULTICOLORGS_HPP

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <functional>
#include <numeric>
#include <tbb/parallel_for.h>

namespace mcl {


/// @brief Multi-Color Gauss-Seidel with a projection operator.
///
/// Implemented for our paper DOI=10.1109/TVCG.2017.2730875
///
/// Examples:
///
/// 
template<typename MatrixType>
class MultiColorGaussSeidel
{
  public:
    using Scalar = typename MatrixType::Scalar;
    using RowSparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::RowMajor>;

    /// @brief Solver options
    struct Options
    {
        int max_iters = 200; ///< max solver iters
        int check_tol = 10;  ///< num iters to check tol
        Scalar rel_tol = Scalar(1e-8);  ///< delta x tolerance
        Scalar omega = Scalar(1);    ///< relaxation parameter
    } options;

    /// @brief Optional: Returns true if the solver should exit, default uses
    /// ||B-AX||/||B|| < options.rel_tol. Example:
    /// is_converged = converged(A, B, X);
    std::function<bool(const RowSparseMatrix&, const MatrixType&, const MatrixType&)> converged;

    /// @brief Optional: Post-sweep projection (e.g., constraints). Example:
    /// project = [&](int i, MatrixXd &X){ X(i,1) = std::max(X(i,1), 0); }; // floor
    /// Called in parallel if colors are used.
    std::function<void(int, MatrixType&)> project;

    /// @brief Solve AX=B s.t. X in C
    /// @param A sparse matrix
    /// @param B linear system rhs
    /// @param X linear system variable
    /// @param colors vertex colors for parallel execution
    /// @param parallel_exec optional flag for parallel/serial color
    /// @return number of iterations taken by solver
    int solve(const RowSparseMatrix& A,
        const MatrixType& B,
        MatrixType& X,
        const std::vector<std::vector<int>>& colors = {},
        const std::vector<bool>& parallel_exec = {});

    /// @brief Checks: ||b-Ax||/||b|| < rel_tol
    bool default_converged(const RowSparseMatrix& A, const MatrixType& B, const MatrixType& X);

    /// @brief Performs an X update at specified row
    void sweep(int row, const RowSparseMatrix& A, const MatrixType& B, MatrixType& X);
};

//
// Implementation
//

    template<typename MatrixType>
int MultiColorGaussSeidel<MatrixType>::solve(const RowSparseMatrix& A,
        const MatrixType& B,
        MatrixType& X,
        const std::vector<std::vector<int>>& colors,
        const std::vector<bool>& parallel_exec)
    {
        int num_colors = colors.size();
        int num_parallel_exec = parallel_exec.size();

        // Resize X to correct size if needed
        if (X.rows() != B.rows() || X.cols() != B.cols()) {
            X = B;
            X.setZero();
        }

        // Set convergence check if not user-provided
        if (converged == nullptr)
        {
            using namespace std::placeholders;
            converged = std::bind(&MultiColorGaussSeidel::default_converged, this, _1, _2, _3);
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
                        sweep(ind, A, B, X);
                        //if (project != nullptr) {
                        //project(ind, X);
                        //}
                    }
                } else {
                    tbb::parallel_for(tbb::blocked_range<int>(0, int(inds.size())),
                                      [&](const tbb::blocked_range<int>& range) {
                                          for (int i = range.begin(); i != range.end(); ++i) {
                                              sweep(inds[i], A, B, X);
                                            if (project != nullptr) {
                                                project(inds[i], X);
                                            }
                                          }
                                      });
                }
            }

            // Check if converged
            if (iter % options.check_tol == 0 && converged(A, B, X)) {
                break;
            }
        }
        return iter;
    }

        template<typename MatrixType>
bool MultiColorGaussSeidel<MatrixType>::default_converged(const RowSparseMatrix& A, const MatrixType& B, const MatrixType& X)
    {
        Scalar B_norm = B.norm();
        Scalar residual = (B - A * X).norm();
        return (residual / B_norm < options.rel_tol);
    }

    template<typename MatrixType>
void MultiColorGaussSeidel<MatrixType>::sweep(int row, const RowSparseMatrix& A, const MatrixType& B, MatrixType& X)
    {
        Scalar Aii = Scalar(0);
        auto LUx = X.row(row); // is there a better way to derive row type?
        LUx.setZero();

        for (typename RowSparseMatrix::InnerIterator iter(A, row); iter; ++iter) {
            if (iter.col() == row) {
                Aii = iter.value();
            } else {
                LUx += iter.value() * X.row(iter.col());
            }
        }
        if (std::abs(Aii) > Scalar(0)) {
            auto X0 = X.row(row);
            auto X1 = (B.row(row) - LUx) / Aii;
            X.row(row) = (Scalar(1) - options.omega) * X0 + options.omega * X1;
        }
    }

} // end namespace mcl

#endif
