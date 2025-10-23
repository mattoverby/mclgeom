// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_XUSPLINE_HPP
#define MCL_GEOM_XUSPLINE_HPP 1

#include <cmath>

namespace mcl {

//
// 	Nonlinear Material Design Using Principal Stretches (2015)
//	Hongyi Xu, Funshing Sin, Yufeng Zhu, Jernej Barbic
//
// The spline class can be overloaded for custom splines.
// Some common ones (StVK, NeoHookean, Co-Rotated linear) are implemented below.
// Code adapted from George Brown (https://georgbrown.github.io)
//
template<typename T>
class XuSpline
{
  public:
    virtual T f(T x) const = 0;
    virtual T g(T x) const = 0;
    virtual T h(T x) const = 0;
    virtual T df(T x) const = 0;
    virtual T dg(T x) const = 0;
    virtual T dh(T x) const = 0;
    virtual T ddf(T x) const = 0;
    virtual T ddg(T x) const = 0;
    virtual T ddh(T x) const = 0;

    // Eq. 16: compression term that helps with stability
    static T compress_term(T kappa, T x) { return (kappa / 12.0) * std::pow((T(1) - x) / 6.0, 3.0); }
    static T d_compress_term(T kappa, T x) { return (-kappa / 24.0) * (std::pow((T(1) - x) / (6.0), 2.0)); }
    static T dd_compress_term(T kappa, T x) { return (T(1) / 432.0) * kappa * (T(1) - x); }
};

template<typename T>
class XuNeoHookean : public XuSpline<T>
{
  public:
    XuNeoHookean(T mu_, T lambda_, T kappa_)
        : mu(mu_)
        , lambda(lambda_)
        , kappa(kappa_)
    {
    }
    const T mu, lambda, kappa;
    T f(T x) const { return 0.5 * mu * (x * x - T(1)); }
    T g(T x) const
    {
        (void)(x);
        return 0.0;
    }
    T h(T x) const
    {
        T logx = std::log(x);
        return -mu * logx + 0.5 * lambda * logx * logx + XuSpline<T>::compress_term(kappa, x);
    }
    T df(T x) const { return mu * x; }
    T dg(T x) const
    {
        (void)(x);
        return 0.0;
    }
    T dh(T x) const { return -mu / x + lambda * std::log(x) / x + XuSpline<T>::d_compress_term(kappa, x); }
    T ddf(T x) const
    {
        (void)(x);
        return mu;
    }
    T ddg(T x) const
    {
        (void)(x);
        return 0.0;
    }
    T ddh(T x) const { return (lambda * (T(1) - std::log(x)) + mu) / (x * x) + XuSpline<T>::dd_compress_term(kappa, x); }
};

template<typename T>
class XuStVK : public XuSpline<T>
{
  public:
    XuStVK(T mu_, T lambda_, T kappa_)
        : mu(mu_)
        , lambda(lambda_)
        , kappa(kappa_)
    {
    }
    const T mu, lambda, kappa;
    T f(T x) const
    {
        T x2 = x * x;
        return 0.125 * lambda * (x2 * x2 - 6.0 * x2 + 5.0) + 0.25 * mu * (x2 - T(1)) * (x2 - T(1));
    }
    T g(T x) const { return 0.25 * lambda * (x * x - T(1)); }
    T h(T x) const { return XuSpline<T>::compress_term(kappa, x); }
    T df(T x) const
    {
        T x2 = x * x;
        return 0.125 * lambda * (4.0 * x2 * x - 12.0 * x) + mu * x * (x2 - T(1));
    }
    T dg(T x) const { return 0.5 * lambda * x; }
    T dh(T x) const { return XuSpline<T>::d_compress_term(kappa, x); }
    T ddf(T x) const { return 0.5 * mu * ((6.0 * x * x - 2.0) + 1.5 * lambda * (x * x - T(1))); }
    T ddg(T x) const
    {
        (void)(x);
        return 0.5 * lambda;
    }
    T ddh(T x) const
    {
        (void)(x);
        return 0.0;
    }
};

template<typename T>
class XuCoRotated : public XuSpline<T>
{
  public:
    XuCoRotated(T mu_, T lambda_, T kappa_)
        : mu(mu_)
        , lambda(lambda_)
        , kappa(kappa_)
    {
    }
    const T mu, lambda, kappa;
    T f(T x) const { return 0.5 * lambda * (x * x - 6.0 * x + 5.0) + mu * (x - T(1)) * (x - T(1)); }
    T g(T x) const { return lambda * (x - T(1)); }
    T h(T x) const { return XuSpline<T>::compress_term(kappa, x); }
    T df(T x) const { return 0.5 * lambda * (2.0 * x - 6.0) + 2.0 * mu * (x - T(1)); }
    T dg(T x) const
    {
        (void)(x);
        return lambda;
    }
    T dh(T x) const { return XuSpline<T>::d_compress_term(kappa, x); }
    T ddf(T x) const
    {
        (void)(x);
        return 2.0 * mu + lambda;
    }
    T ddg(T x) const
    {
        (void)(x);
        return 0.0;
    }
    T ddh(T x) const
    {
        (void)(x);
        return 0.0;
    }
};

} // end namespace mcl

#endif
