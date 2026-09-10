// Copyright Matt Overby 2026.
// Distributed under the MIT License.

#ifndef MCL_GEOM_LMSOLVER_HPP
#define MCL_GEOM_LMSOLVER_HPP

#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <functional>

namespace mcl {

/// @brief Levenberg-Marquard solver for underdetermined systems.
/// Follows derivations from https://doi.org/10.1111/cgf.14361, Eq. 12.
/// Solves (1/2)||f(x)||^2, with f : R^n -> R^m
/// Each iteration computes step: p = -J^T(JJ^T + U)*f(x)
/// for Jacobian J, damping matrix U, and residual f(x).
template<typename VectorType, typename SparseMatrixType>
class LevenbergMarquardt
{
  public:
    using T = typename VectorType::Scalar;
    using LDLT = Eigen::SimplicialLDLT<SparseMatrixType>;

    struct Options
    {
        int max_iters = 30;       ///< max solver iters
        T tol = T(1e-6);          ///< convergence tol
        T lm_param = T(1e-4);     ///< diagonal regularizer, adjusted in iterate(x)
        T min_lm_param = T(1e-8); ///< if nonnegative, enables adaptive LM
        int max_ls_iters = 1000;  ///< linesearch iterations
    } options;

    /// @brief Required: computes residual and Jacobian J = grad f(x)
    /// The Jacobian isn't needed if the last argument is false (i.e., linesearch).
    /// Signature is [x, residual, Jacobian, needsJacobian]
    std::function<void(const VectorType&, VectorType&, SparseMatrixType&, bool)> objective;

    /// @brief Performs an LM iteration of f(x), f : R^n -> R^m
    /// Returns the objective (1/2)||f(x)||^2 and updates x, or -1 if there was an error.
    T iterate(VectorType& x)
    {
        // Compute Jacobian
        SparseMatrixType J;
        VectorType residual;
        objective(x, residual, J, true);
        if (J.nonZeros() == 0) {
            return T(0);
        }

        VectorType lambda = VectorType::Zero(J.rows());
        auto JJt = (J * J.transpose()).eval();
        T eval_init = T(0.5) * residual.dot(residual);

        // Special cases: small J use dense solvers
        // Solve y = JJ^T(f(x))
        if (J.rows() > 3) {
            T max_JJti = JJt.diagonal().maxCoeff();
            for (int i = 0; i < JJt.rows(); ++i) {
                JJt.coeffRef(i, i) += max_JJti * options.lm_param;
            }
            // JJt.diagonal().array() += max_JJti * options.lm_param;
            LDLT ldlt(JJt);
            if (ldlt.info() == Eigen::Success) {
                lambda = ldlt.solve(residual);
            }
        } else if (J.rows() == 1) {
            T JJt0 = JJt.coeff(0, 0);
            lambda[0] = residual[0] / JJt0;
        } else if (J.rows() == 2) {
            Eigen::Matrix<T, 2, 2> JJt_dense(JJt);
            Eigen::JacobiSVD<Eigen::Matrix<T, 2, 2>> svd(JJt_dense, Eigen::ComputeFullU | Eigen::ComputeFullV);
            lambda = svd.solve(residual);
        } else if (J.rows() == 3) {
            Eigen::Matrix<T, 3, 3> JJt_dense(JJt);
            Eigen::JacobiSVD<Eigen::Matrix<T, 3, 3>> svd(JJt_dense, Eigen::ComputeFullU | Eigen::ComputeFullV);
            lambda = svd.solve(residual);
        }

        // Search direction
        VectorType p = J.transpose() * (-lambda);
        T p_norm = p.template lpNorm<Eigen::Infinity>();
        if (!std::isfinite(p_norm)) {
            return T(-1);
        }

        SparseMatrixType dummy;
        T alpha = T(1); // step size

        // Initial step
        VectorType x0 = x;
        x = x0 + alpha * p;
        objective(x, residual, dummy, false);
        T eval_new = T(0.5) * residual.dot(residual);

        // Linesearch
        int ls_iter = 0;
        while (eval_new > eval_init) {
            alpha *= T(0.5);
            x = x0 + alpha * p;
            objective(x, residual, dummy, false);
            eval_new = T(0.5) * residual.dot(residual);
            ls_iter++;
            if (ls_iter > options.max_ls_iters) {
                eval_new = T(-1); // did not converge
                break;
            }
        }

        // Update LM damping parameter
        if (eval_new >= 0 && options.min_lm_param >= 0) {
            if (alpha < 0.25) {
                options.lm_param *= T(2);
            } // take a smaller step
            else if (alpha > 0.75) {
                options.lm_param /= T(3);
            } // take a bigger step
            options.lm_param = std::clamp(options.lm_param, options.min_lm_param, T(1e-1));
        }

        return eval_new;
    }
};

} // end namespace mcl

#endif // MCL_GEOM_LMSOLVER_HPP
