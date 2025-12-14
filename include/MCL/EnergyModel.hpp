// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_ENERGYMODEL_HPP
#define MCL_GEOM_ENERGYMODEL_HPP

#include "Lame.hpp"
#include "XuSpline.hpp"
#include <Eigen/Dense>

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
// See: MCL/XuSpline.hpp
template<int DIM, typename T>
class XuSplineModel
{
  protected:
    typedef Eigen::Matrix<T, DIM, 1> VecD;
    typedef Eigen::Matrix<T, DIM, DIM> MatD;

  public:
    static T energy_density(const XuSpline<T>* s, const VecD& x);
    static T gradient(const XuSpline<T>* s, const VecD& x, VecD& g);
    static void hessian(const XuSpline<T>* s, const VecD& x, MatD& H);
};

// Simplified Stable Neo-Hookean (similar to Smith et al. '17)
template<int DIM, typename T>
class StableNeoHookean
{
  protected:
    typedef Eigen::Matrix<T, DIM, 1> VecD;
    typedef Eigen::Matrix<T, DIM, DIM> MatD;

  public:
    static T energy_density(const Lame<T>& lame, const VecD& x);
    static T gradient(const Lame<T>& lame, const VecD& x, VecD& g);
    static void hessian(const Lame<T>& lame, const VecD& x, MatD& H);
};

// Classic (isotropic) Neo-Hookean
template<int DIM, typename T>
class IsoNeoHookean
{
  protected:
    typedef Eigen::Matrix<T, DIM, 1> VecD;
    typedef Eigen::Matrix<T, DIM, DIM> MatD;

  public:
    static T energy_density(const Lame<T>& lame, const VecD& x);
    static T gradient(const Lame<T>& lame, const VecD& x, VecD& g);
    static void hessian(const Lame<T>& lame, const VecD& x, MatD& H);
};

// SLIM Symmetric Dirichlet
// f(x) = (||F||^2 + ||F^-1||^2)/2
template<int DIM, typename T>
class SymmDirichlet
{
  protected:
    typedef Eigen::Matrix<T, DIM, 1> VecD;
    typedef Eigen::Matrix<T, DIM, DIM> MatD;

  public:
    static T energy_density(const VecD& x);
    static T gradient(const VecD& x, VecD& g);
    static void hessian(const VecD& x, MatD& H);
};

//
//	Implementation
//

//================================================
//	XuSpline
//================================================

template<int DIM, typename T>
T
XuSplineModel<DIM, T>::energy_density(const XuSpline<T>* s, const VecD& x)
{
    if (DIM == 3) {
        return s->f(x[0]) + s->f(x[1]) + s->f(x[2]) + s->g(x[0] * x[1]) + s->g(x[1] * x[2]) + s->g(x[2] * x[0]) +
               s->h(x[0] * x[1] * x[2]);
    }
    return s->f(x[0]) + s->f(x[1]) + s->g(x[0] * x[1]) + s->h(x[0] * x[1]);
}

template<int DIM, typename T>
T
XuSplineModel<DIM, T>::gradient(const XuSpline<T>* s, const VecD& x, VecD& g)
{
    if (DIM == 3) {
        T hprime = s->dh(x[0] * x[1] * x[2]);
        g[0] = s->df(x[0]) + s->dg(x[0] * x[1]) * x[1] + s->dg(x[2] * x[0]) * x[2] + hprime * x[1] * x[2];
        g[1] = s->df(x[1]) + s->dg(x[1] * x[2]) * x[2] + s->dg(x[0] * x[1]) * x[0] + hprime * x[2] * x[0];
        g[2] = s->df(x[2]) + s->dg(x[2] * x[0]) * x[0] + s->dg(x[1] * x[2]) * x[1] + hprime * x[0] * x[1];
        return energy_density(s, x);
    }
    T hprime = s->dh(x[0] * x[1]);
    g[0] = s->df(x[0]) + s->dg(x[0] * x[1]) * x[1] + hprime * x[1];
    g[1] = s->df(x[1]) + s->dg(x[0] * x[1]) * x[0] + hprime * x[0];
    return energy_density(s, x);
}

template<int DIM, typename T>
void
XuSplineModel<DIM, T>::hessian(const XuSpline<T>* s, const VecD& x, MatD& H)
{
    if (DIM == 3) {
        T hprime = s->dh(x[0] * x[1] * x[2]);
        T hprimeprime = s->ddh(x[0] * x[1] * x[2]);
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
    T hprimeprime = s->ddh(x[0] * x[1]);
    H(0, 0) = s->ddf(x[0]) + s->ddg(x[0] * x[1]) * x[1] * x[1] + hprimeprime * x[1] * x[1];
    H(1, 1) = s->ddf(x[1]) + s->ddg(x[0] * x[1]) * x[0] * x[0] + hprimeprime * x[0] * x[0];
    H(0, 1) = s->ddg(x[0] * x[1]) * x[0] * x[1] + s->dg(x[0] * x[1]) + hprimeprime * x[0] * x[1];
    H(1, 0) = H(0, 1);
}

//================================================
//	Stable Neo-Hookean
//================================================

template<int DIM, typename T>
T
StableNeoHookean<DIM, T>::energy_density(const Lame<T>& lame, const VecD& x)
{
    if (!x.allFinite() || x.minCoeff() <= 0)
        return std::numeric_limits<float>::max();

    T alpha = 1.0 + lame.mu() / lame.lambda();
    constexpr T d = T(DIM);
    T I_1 = x.squaredNorm();
    T J = x.prod();
    return 0.5 * lame.lambda() * (J - alpha) * (J - alpha) + 0.5 * lame.mu() * (I_1 - d);
}

template<int DIM, typename T>
T
StableNeoHookean<DIM, T>::gradient(const Lame<T>& lame, const VecD& x, VecD& grad)
{
    T alpha = 1.0 + lame.mu() / lame.lambda();
    T J = x.prod();
    VecD x_inv = x.cwiseInverse();
    grad = lame.mu() * x + lame.lambda() * (J - alpha) * J * x_inv;
    return energy_density(lame, x);
}

template<int DIM, typename T>
void
StableNeoHookean<DIM, T>::hessian(const Lame<T>& lame, const VecD& x, MatD& hess)
{
    T alpha = 1.0 + lame.mu() / lame.lambda();
    T J = x.prod();
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

template<int DIM, typename T>
T
IsoNeoHookean<DIM, T>::energy_density(const Lame<T>& lame, const VecD& x)
{
    if (!x.allFinite() || x.minCoeff() <= 0)
        return std::numeric_limits<float>::max();

    T J = x.prod();
    T I_1 = x.squaredNorm();
    T I_3 = J * J;
    T log_I3 = std::log(I_3);
    T t1 = 0.5 * lame.mu() * (I_1 - log_I3 - 3.0);
    T t2 = 0.125 * lame.lambda() * log_I3 * log_I3;
    T r = t1 + t2;
    return r;
}

template<int DIM, typename T>
T
IsoNeoHookean<DIM, T>::gradient(const Lame<T>& lame, const VecD& x, VecD& grad)
{
    T J = x.prod();
    VecD x_inv = x.cwiseInverse();
    grad = (lame.mu() * (x - x_inv) + lame.lambda() * std::log(J) * x_inv);
    return energy_density(lame, x);
}

template<int DIM, typename T>
void
IsoNeoHookean<DIM, T>::hessian(const Lame<T>& lame, const VecD& x, MatD& hess)
{
    static const MatD Iden = MatD::Identity();
    T J = x.prod();
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

template<int DIM, typename T>
T
SymmDirichlet<DIM, T>::energy_density(const VecD& x)
{
    T e = 0;
    for (int i = 0; i < DIM; ++i)
        e += x[i] * x[i] + std::pow(x[i], T(-2.0));

    return e;
}

template<int DIM, typename T>
T
SymmDirichlet<DIM, T>::gradient(const VecD& x, VecD& grad)
{
    T e = 0;
    VecD pow_x;
    for (int i = 0; i < DIM; ++i) {
        pow_x[i] = std::pow(x[i], T(-3.0));
        e += x[i] * x[i] + std::pow(x[i], T(-2.0));
    }
    grad = 2.0 * (x - pow_x);
    return e;
}

template<int DIM, typename T>
void
SymmDirichlet<DIM, T>::hessian(const VecD& x, MatD& H)
{
    H.setZero();
    for (int i = 0; i < DIM; ++i)
        H(i, i) = 2.0 + 6.0 * std::pow(x[i], T(-4.0));
}

//================================================
// Wrapper/catch-all for known energy models.
//================================================

template<int DIM, typename T>
class EnergyModel
{
  protected:
    typedef Eigen::Matrix<T, DIM, 1> VecD;
    typedef Eigen::Matrix<T, DIM, DIM> MatD;

  public:
    // Should consider a more efficient way to handle spline model.
    // It's an object instead of a static class to allow custom splines to
    // be used XuSplineModel. So for now declare known ones as objects.
    const Lame<T>& lame;
    const XuNeoHookean<T> xu_nh;
    const XuStVK<T> xu_stvk;
    const XuCoRotated<T> xu_cr;

    EnergyModel(const Lame<T>& lame_)
        : lame(lame_)
        , xu_nh(lame_.mu(), lame_.lambda(), 0.1)
        , xu_stvk(lame_.mu(), lame_.lambda(), 0.1)
        , xu_cr(lame_.mu(), lame_.lambda(), 0.1)
    {
    }

    T energy_density(const VecD& x)
    {
        T e = std::numeric_limits<T>::max();
        switch (lame.model()) {
            case ENERGY_MODEL_STABLE_NH: {
                e = StableNeoHookean<DIM, T>::energy_density(lame, x);
            } break;
            case ENERGY_MODEL_SYMM_DIR: {
                e = SymmDirichlet<DIM, T>::energy_density(x);
            } break;
            case ENERGY_MODEL_NH: {
                e = IsoNeoHookean<DIM, T>::energy_density(lame, x);
            } break;
            case ENERGY_MODEL_XUSPLINE_NH: {
                e = XuSplineModel<DIM, T>::energy_density(&xu_nh, x);
            } break;
            case ENERGY_MODEL_XUSPLINE_STVK: {
                e = XuSplineModel<DIM, T>::energy_density(&xu_stvk, x);
            } break;
            case ENERGY_MODEL_XUSPLINE_COROTATE: {
                e = XuSplineModel<DIM, T>::energy_density(&xu_cr, x);
            } break;
            default:
                break;
        }
        return e;
    }

    T gradient(const VecD& x, VecD& g)
    {
        T e = std::numeric_limits<T>::max();
        switch (lame.model()) {
            case ENERGY_MODEL_STABLE_NH: {
                e = StableNeoHookean<DIM, T>::gradient(lame, x, g);
            } break;
            case ENERGY_MODEL_SYMM_DIR: {
                e = SymmDirichlet<DIM, T>::gradient(x, g);
            } break;
            case ENERGY_MODEL_NH: {
                e = IsoNeoHookean<DIM, T>::gradient(lame, x, g);
            } break;
            case ENERGY_MODEL_XUSPLINE_NH: {
                e = XuSplineModel<DIM, T>::gradient(&xu_nh, x, g);
            } break;
            case ENERGY_MODEL_XUSPLINE_STVK: {
                e = XuSplineModel<DIM, T>::gradient(&xu_stvk, x, g);
            } break;
            case ENERGY_MODEL_XUSPLINE_COROTATE: {
                e = XuSplineModel<DIM, T>::gradient(&xu_cr, x, g);
            } break;
            default:
                break;
        }
        return e;
    }

    void hessian(const VecD& x, MatD& H)
    {
        H.setZero();
        switch (lame.model()) {
            case ENERGY_MODEL_STABLE_NH: {
                StableNeoHookean<DIM, T>::hessian(lame, x, H);
            } break;
            case ENERGY_MODEL_SYMM_DIR: {
                SymmDirichlet<DIM, T>::hessian(x, H);
            } break;
            case ENERGY_MODEL_NH: {
                IsoNeoHookean<DIM, T>::hessian(lame, x, H);
            } break;
            case ENERGY_MODEL_XUSPLINE_NH: {
                XuSplineModel<DIM, T>::hessian(&xu_nh, x, H);
            } break;
            case ENERGY_MODEL_XUSPLINE_STVK: {
                XuSplineModel<DIM, T>::hessian(&xu_stvk, x, H);
            } break;
            case ENERGY_MODEL_XUSPLINE_COROTATE: {
                XuSplineModel<DIM, T>::hessian(&xu_cr, x, H);
            } break;
            default:
                break;
        }
    }
};

} // end namespace mcl

#endif
