// Copyright Matt Overby 2025.
// Distributed under the MIT License.

#ifndef MCL_GEOM_KKTSOLVER_HPP
#define MCL_GEOM_KKTSOLVER_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <functional>

namespace mcl {

/// @brief Conjugate Gradient for saddle point system (KKT matrix):
/// [ A  C^T ] [ x ] = [ b ]
/// [ C  0   ] [ y ] = [ d ]
/// From DOI 10.1109/TVCG.2017.2730875, Alg. 2
template<typename VectorType, typename SparseMatrixType>
class KKTSolver
{
  protected:
    using Scalar = typename VectorType::Scalar;
    using LDLT = Eigen::SimplicialLDLT<SparseMatrixType>;
    VectorType q1, q2, q3, r, s;

  public:
    struct Options
    {
        int max_iters = 30; ///< max solver iters
        Scalar tol = 1e-6;  ///< convergence tol
    } options;

    /// @brief Optional: computes x = A^-1(b).
    std::function<void(const VectorType&, VectorType&)> solve_Axb;

    /// @brief Resizes internal cache variables used during solve.
    void resize(int A_rows, int C_rows, int C_cols)
    {
        q1.resize(A_rows);
        q2.resize(A_rows);
        q3.resize(C_rows);
        r.resize(C_rows);
        s.resize(C_rows);
        q1.setZero();
        q2.setZero();
        q3.setZero();
        r.setZero();
        s.setZero();
    }

    /// @brief Solve for x and y, return number of iterations.
    /// A is allowed to be empty if solve_Axb function is defined.
    /// Otherwise, it is factorized with Cholesky.
    int solve(const SparseMatrixType& A,
              const VectorType& b,
              const SparseMatrixType& C,
              const VectorType& d,
              VectorType& x,
              VectorType& y)
    {
        // Factorize A if we don't have a function to get x.
        auto ldlt_ptr = std::make_shared<LDLT>();
        if (solve_Axb == nullptr) {
            ldlt_ptr->compute(A);
            if (ldlt_ptr->info() != Eigen::Success) {
                return -1;
            }
            solve_Axb = [ldlt_ptr](const VectorType& bt, VectorType& xt) { xt = ldlt_ptr->solve(bt); };
        }

        // No constraints, just use LDLT
        if (C.nonZeros() == 0) {
            solve_Axb(b, x);
            return 1;
        }

        resize(A.rows(), C.rows(), C.cols());

        if (y.rows() != C.rows()) {
            y.resize(C.rows());
            y.setZero();
        }

        auto Ct = C.transpose();
        q1 = b - Ct * y;
        solve_Axb(q1, x);

        r = C * x - d;
        s = r.eval();
        q2.setZero();
        q3.setZero();
        Scalar tol2 = options.tol * options.tol;

        int iter = 0;
        while (iter < options.max_iters) {
            iter++;

            // Intermediate variables of Schur compliment:
            q1 = Ct * s;
            solve_Axb(q1, q2);
            q3 = C * q2;

            Scalar denom = s.dot(q3);
            if (denom == Scalar(0)) {
                break;
            }
            Scalar alpha = s.dot(r) / denom;

            // Take a step
            x -= alpha * q2;
            y += alpha * s;
            r -= alpha * q3;

            // Exit if resid is low enough
            if (r.squaredNorm() < tol2) {
                break;
            }

            denom = s.dot(q3);
            if (denom == Scalar(0)) {
                break;
            }
            Scalar beta = r.dot(q3) / denom;
            s = r - (beta * s);
        }

        return iter;
    }
};

} // end namespace mcl

#endif // MCL_GEOM_KKTSOLVER_HPP
