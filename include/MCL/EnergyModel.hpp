// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_ENERGYMODEL_HPP
#define MCL_ENERGYMODEL_HPP 1

#include "MCL/Lame.hpp"
#include "MCL/XuSpline.hpp"
#include <Eigen/Dense>
#include <memory>

namespace mcl {

enum
{
    ENERGY_MODEL_INVALID = -1,
    ENERGY_MODEL_ARAP,
    ENERGY_MODEL_XUSPLINE_NH,       // Xu spline (neo-Hookean)
    ENERGY_MODEL_XUSPLINE_STVK,     // Xu spline (St. VK)
    ENERGY_MODEL_XUSPLINE_COROTATE, // Xu spline (CR)
    ENERGY_MODEL_STABLE_NH,         // Stable neo-Hookean
    ENERGY_MODEL_NH,                // iso log-barrier neo-Hookean
    ENERGY_MODEL_SYMM_DIR,          // Symmetric Dirichlet
    ENERGY_MODEL_NUM
};

// Nonlinear Material Design Using Principal Stretches, Xu et al. 2015.
template<int DIM>
class XuSplineModel
{
  protected:
    typedef Eigen::Matrix<double, DIM, 1> VecD;
    typedef Eigen::Matrix<double, DIM, DIM> MatD;

  public:
    static double energy_density(const XuSpline<double>* s, const VecD& x);
    static double gradient(const XuSpline<double>* s, const VecD& x, VecD& g);
    static void hessian(const XuSpline<double>* s, const VecD& x, MatD& H);
};

// Simplified stable Neo-Hookean (similar to Smith et al. '17)
template<int DIM>
class StableNeoHookean
{
  protected:
    typedef Eigen::Matrix<double, DIM, 1> VecD;
    typedef Eigen::Matrix<double, DIM, DIM> MatD;

  public:
    static double energy_density(const Lame& lame, const VecD& x);
    static double gradient(const Lame& lame, const VecD& x, VecD& g);
    static void hessian(const Lame& lame, const VecD& x, MatD& H);
};

// Classic (isotropic) Neo-Hookean
template<int DIM>
class IsoNeoHookean
{
  protected:
    typedef Eigen::Matrix<double, DIM, 1> VecD;
    typedef Eigen::Matrix<double, DIM, DIM> MatD;

  public:
    static double energy_density(const Lame& lame, const VecD& x);
    static double gradient(const Lame& lame, const VecD& x, VecD& g);
    static void hessian(const Lame& lame, const VecD& x, MatD& H);
};

// SLIM Symmetric Dirichlet
// f(x) = (||F||^2 + ||F^-1||^2)/2
template<int DIM>
class SymmDirichlet
{
  protected:
    typedef Eigen::Matrix<double, DIM, 1> VecD;
    typedef Eigen::Matrix<double, DIM, DIM> MatD;

  public:
    static double energy_density(const VecD& x);
    static double gradient(const VecD& x, VecD& g);
    static void hessian(const VecD& x, MatD& H);
};

//
//	Implementation
//

//================================================
//	XuSpline
//================================================

template<int DIM>
double
XuSplineModel<DIM>::energy_density(const XuSpline<double>* s, const VecD& x)
{
    if (DIM == 3) {
        return s->f(x[0]) + s->f(x[1]) + s->f(x[2]) + s->g(x[0] * x[1]) + s->g(x[1] * x[2]) + s->g(x[2] * x[0]) +
               s->h(x[0] * x[1] * x[2]);
    }
    return s->f(x[0]) + s->f(x[1]) + s->g(x[0] * x[1]) + s->h(x[0] * x[1]);
}

template<int DIM>
double
XuSplineModel<DIM>::gradient(const XuSpline<double>* s, const VecD& x, VecD& g)
{
    if (g.size() != x.size()) {
        g.resize(x.size());
    }
    if (DIM == 3) {
        double hprime = s->dh(x[0] * x[1] * x[2]);
        g[0] = s->df(x[0]) + s->dg(x[0] * x[1]) * x[1] + s->dg(x[2] * x[0]) * x[2] + hprime * x[1] * x[2];
        g[1] = s->df(x[1]) + s->dg(x[1] * x[2]) * x[2] + s->dg(x[0] * x[1]) * x[0] + hprime * x[2] * x[0];
        g[2] = s->df(x[2]) + s->dg(x[2] * x[0]) * x[0] + s->dg(x[1] * x[2]) * x[1] + hprime * x[0] * x[1];
        return energy_density(s, x);
    }
    double hprime = s->dh(x[0] * x[1]);
    g[0] = s->df(x[0]) + s->dg(x[0] * x[1]) * x[1] + hprime * x[1];
    g[1] = s->df(x[1]) + s->dg(x[0] * x[1]) * x[0] + hprime * x[0];
    return energy_density(s, x);
}

template<int DIM>
void
XuSplineModel<DIM>::hessian(const XuSpline<double>* s, const VecD& x, MatD& H)
{
    if (DIM == 3) {
        double hprime = s->dh(x[0] * x[1] * x[2]);
        double hprimeprime = s->ddh(x[0] * x[1] * x[2]);
        H(0, 0) = s->ddf(x[0]) + s->ddg(x[0] * x[1]) * x[1] * x[1] + s->ddg(x[2] * x[0]) * x[2] * x[2] +
                  hprimeprime * x[1] * x[1] * x[2] * x[2];
        H(1, 1) = s->ddf(x[1]) + s->ddg(x[1] * x[2]) * x[2] * x[2] + s->ddg(x[0] * x[1]) * x[0] * x[0] +
                  hprimeprime * x[2] * x[2] * x[0] * x[0];
        H(2, 2) = s->ddf(x[2]) + s->ddg(x[2] * x[0]) * x[0] * x[0] + s->ddg(x[1] * x[2]) * x[1] * x[1] +
                  hprimeprime * x[0] * x[0] * x[1] * x[1];
        H(0, 1) = s->ddg(x[0] * x[1]) * x[0] * x[1] + s->dg(x[0] * x[1]) + hprimeprime * x[0] * x[1] * x[2] * x[2] +
                  hprime * x[2];
        H(1, 0) = H(0, 1);
        H(0, 2) = s->ddg(x[0] * x[2]) * x[0] * x[2] + s->dg(x[0] * x[2]) + hprimeprime * x[0] * x[2] * x[1] * x[1] +
                  hprime * x[1];
        H(2, 0) = H(0, 2);
        H(1, 2) = s->ddg(x[1] * x[2]) * x[1] * x[2] + s->dg(x[1] * x[2]) + hprimeprime * x[1] * x[2] * x[0] * x[0] +
                  hprime * x[0];
        H(2, 1) = H(1, 2);
        return;
    }
    double hprimeprime = s->ddh(x[0] * x[1]);
    H(0, 0) = s->ddf(x[0]) + s->ddg(x[0] * x[1]) * x[1] * x[1] + hprimeprime * x[1] * x[1];
    H(1, 1) = s->ddf(x[1]) + s->ddg(x[0] * x[1]) * x[0] * x[0] + hprimeprime * x[0] * x[0];
    H(0, 1) = s->ddg(x[0] * x[1]) * x[0] * x[1] + s->dg(x[0] * x[1]) + hprimeprime * x[0] * x[1];
    H(1, 0) = H(0, 1);
}

//================================================
//	Stable Neo-Hookean
//================================================

template<int DIM>
double
StableNeoHookean<DIM>::energy_density(const Lame& lame, const VecD& x)
{
    if (!x.allFinite() || x.minCoeff() <= 0)
        return std::numeric_limits<float>::max();

    double alpha = 1.0 + lame.mu() / lame.lambda();
    constexpr double d = double(DIM);
    double I_1 = x.squaredNorm();
    double J = x.prod();
    return 0.5 * lame.lambda() * (J - alpha) * (J - alpha) + 0.5 * lame.mu() * (I_1 - d);
}

template<int DIM>
double
StableNeoHookean<DIM>::gradient(const Lame& lame, const VecD& x, VecD& grad)
{
    double alpha = 1.0 + lame.mu() / lame.lambda();
    double J = x.prod();
    VecD x_inv = x.cwiseInverse();
    grad = lame.mu() * x + lame.lambda() * (J - alpha) * J * x_inv;
    return energy_density(lame, x);
}

template<int DIM>
void
StableNeoHookean<DIM>::hessian(const Lame& lame, const VecD& x, MatD& hess)
{
    double alpha = 1.0 + lame.mu() / lame.lambda();
    double J = x.prod();
    if (DIM == 3) {
        VecD coeffs;
        coeffs[0] = x[1] * x[2];
        coeffs[1] = x[0] * x[2];
        coeffs[2] = x[0] * x[1];
        MatD coeffsMat;
        coeffsMat.setZero();
        coeffsMat(0, 1) = x[2];
        coeffsMat(1, 0) = x[2];
        coeffsMat(0, 2) = x[1];
        coeffsMat(2, 0) = x[1];
        coeffsMat(1, 2) = x[0];
        coeffsMat(2, 1) = x[0];
        hess = lame.mu() * MatD::Identity() + lame.lambda() * coeffs * coeffs.transpose() +
               lame.lambda() * (J - alpha) * coeffsMat;
        return;
    }
    hess(0, 0) = lame.lambda() * x[1] * x[1] + lame.mu();
    hess(1, 1) = lame.lambda() * x[0] * x[0] + lame.mu();
    hess(0, 1) = lame.lambda() * (J - alpha) + lame.lambda() * J;
    hess(1, 0) = hess(0, 1);
}

//================================================
//	Isotropic Neo-Hookean
//================================================

template<int DIM>
double
IsoNeoHookean<DIM>::energy_density(const Lame& lame, const VecD& x)
{
    if (!x.allFinite() || x.minCoeff() <= 0)
        return std::numeric_limits<float>::max();

    double J = x.prod();
    double I_1 = x.squaredNorm();
    double I_3 = J * J;
    double log_I3 = std::log(I_3);
    double t1 = 0.5 * lame.mu() * (I_1 - log_I3 - 3.0);
    double t2 = 0.125 * lame.lambda() * log_I3 * log_I3;
    double r = t1 + t2;
    return r;
}

template<int DIM>
double
IsoNeoHookean<DIM>::gradient(const Lame& lame, const VecD& x, VecD& grad)
{
    double J = x.prod();
    VecD x_inv = x.cwiseInverse();
    grad = (lame.mu() * (x - x_inv) + lame.lambda() * std::log(J) * x_inv);
    return energy_density(lame, x);
}

template<int DIM>
void
IsoNeoHookean<DIM>::hessian(const Lame& lame, const VecD& x, MatD& hess)
{
    static const MatD Iden = MatD::Identity();
    double J = x.prod();
    VecD x_inv = x.cwiseInverse();
    MatD invXmat;
    invXmat.setZero();
    invXmat(0, 0) = 1.0 / (x[0] * x[0]);
    invXmat(1, 1) = 1.0 / (x[1] * x[1]);
    if (DIM == 3) {
        invXmat(2, 2) = 1.0 / (x[2] * x[2]);
    }
    hess = lame.mu() * (Iden - 2.0 * invXmat) + lame.lambda() * std::log(J) * invXmat +
           lame.lambda() * x_inv * x_inv.transpose();
}

//================================================
//	Symmetric Dirichlet
//================================================

template<int DIM>
double
SymmDirichlet<DIM>::energy_density(const VecD& x)
{
    double e = 0;
    for (int i = 0; i < DIM; ++i)
        e += x[i] * x[i] + std::pow(x[i], -2.0);

    return e;
}

template<int DIM>
double
SymmDirichlet<DIM>::gradient(const VecD& x, VecD& grad)
{
    double e = 0;
    VecD pow_x;
    for (int i = 0; i < DIM; ++i) {
        pow_x[i] = std::pow(x[i], -3.0);
        e += x[i] * x[i] + std::pow(x[i], -2.0);
    }
    grad = 2.0 * (x - pow_x);
    return e;
}

template<int DIM>
void
SymmDirichlet<DIM>::hessian(const VecD& x, MatD& H)
{
    H.setZero();
    for (int i = 0; i < DIM; ++i)
        H(i, i) = 2.0 + 6.0 * std::pow(x[i], -4.0);
}

} // end namespace mcl

#endif
