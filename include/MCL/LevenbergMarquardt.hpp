// Copyright Matt Overby 2026.
// Distributed under the MIT License.

#ifndef MCL_GEOM_LMSOLVER_HPP
#define MCL_GEOM_LMSOLVER_HPP

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/JacobiSVD>
#include <functional>

namespace mcl {

/// @brief Levenberg-Marquard linear system solver for underdetermined systems.
/// Follows derivations from https://doi.org/10.1111/cgf.14361, Eq. 12.
/// Solves (1/2)||f(x)||^2, with f : R^n -> R^m
template<typename VectorType, typename SparseMatrixType>
class LevenbergMarquardt
{
  public:
    using Scalar = typename VectorType::Scalar;
    using LDLT = Eigen::SimplicialLDLT<SparseMatrixType>;

    struct Options
    {
        int max_iters = 30; ///< max solver iters
        Scalar tol = 1e-6;  ///< convergence tol
        Scalar lm_param = 1e-4; ///< diagonal regularizer
        Scalar min_lm_param = 1e-8; ///< if nonnegative, enables adaptive LM
        int max_ls_iters = 1000;
    } options;

    /// @brief Required: computes residual and Jacobian J = grad f(x)
    /// The Jacobian isn't needed if the last argument is false (i.e., linesearch).
    std::function<void(const VectorType&, VectorType&, SparseMatrixType&, bool)> objective;

    /// @brief Performs an LM iteration of f(x), f : R^n -> R^m
    /// Returns the objective (1/2)||f(x)||^2 and updates x, or -1 if there was an error.
    T iterate(VectorType &x)
    {
        // Compute Jacobian
        SparseMatrixType J;
        VectorType residual;
        objective(x, J, residual, true);
        if (J.rows() == 0) {
            return T(0);
        }

        VectorType lambda = VectorType::Zero(J.rows());
        auto JJt = J*J.transpose();
        T eval_init = T(0.5) * residual.dot(residual);

        // Special cases: small J use dense solvers
        // Solve y = JJ^T(f(x))
        if (J.rows() > 3) {
            T max_JJti = JJt.diagonal().maxCoeff();
            JJt.diagonal().array() += max_JJti * options.lm_param;
            LDLT ldlt(JJt);
            if(ldlt.info() == Eigen::Success) {
                lambda = ldlt.solve(residual);
            }
        }
        else if (jr == 1) {
            T JJt0 = JJt.coeff(0,0);
            lambda[0] = residual[0] / JJt0;
        }
        else if (jr == 2) {
            Eigen::Matrix<T,2,2> JJt_dense(JJt);
            Eigen::JacobiSVD<Eigen::Matrix<T,2,2>> svd(JJt_dense, ComputeFullU | ComputeFullV);
            lambda = svd.solve(residual);
        }
        else if (jr == 3) {
            Eigen::Matrix<T,3,3> JJt_dense(JJt);
            Eigen::JacobiSVD<Eigen::Matrix<T,3,3>> svd(JJt_dense, ComputeFullU | ComputeFullV);
            lambda = svd.solve(residual);
        }
    
        // Search direction
        VectorType p = J.transpose()*(-lambda);
        T p_norm = p.lpNorm<Infinity>();
        if (!std::isfinite(p_norm)) {
            return T(-1);
        }

        SparseMatrixType dummy;
        T alpha = T(1); // step size

        // Initial step
        VectorType x0 = x;
        x = x0 + alpha * p;
        objective(x, dummy, residual, false);
        T eval_new = T(0.5) * residual.dot(residual);

        // Linesearch
        while (eval_new > eval_init) {
            alpha *= T(0.5);
            x = x0 + alpha * p;
            objective(x, dummy, residual, false);
            eval_new = T(0.5) * residual.dot(residual);
            ls_iter++;
            if(ls_iter > options.max_ls_iters){
                eval_new = T(-1); // did not converge
                break;
            }
        }
    

        // Update LM damping parameter
        if (eval_new >= 0 && options.min_lm_param >= 0) {
            if (alpha < 0.25) { lm_param *= 2.0; } // take a smaller step
            else if (alpha > 0.75) { lm_param /= 3.0; } // take a bigger step
            lm_param = std::clamp(lm_param, options.min_lm_param, 1e-1);
        }
        
        return eval_new;
}

};

} // end namespace mcl

#endif // MCL_GEOM_LMSOLVER_HPP
