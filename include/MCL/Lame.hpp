// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_LAME_HPP
#define MCL_GEOM_LAME_HPP 1

#include <cmath>

namespace mcl {

template <typename T>
class Lame
{
  public:
    // Static presets
    static Lame rubber(int model_ = 0) { return Lame(10000000, 0.499, model_); }          // true rubber
    static Lame soft_rubber(int model_ = 0) { return Lame(10000000, 0.399, model_); }     // fun rubber!
    static Lame very_soft_rubber(int model_ = 0) { return Lame(1000000, 0.299, model_); } // more fun!
    static Lame one(int model_ = 0)
    {
        Lame l;
        l.set_lame(0, 1);
        l.set_model(model_);
        return l;
    } // bulkmod = 1

    // k: Youngs (Pa), measure of stretch
    // v: Poisson, measure of incompressibility
    Lame(T k, T v, int model_ = 0)
        : m_limit_min(-100.0)
        , m_limit_max(100.0)
    {
        set_yp(k, v);
        set_model(model_);
    }

    // Default = soft rubber
    Lame()
        : Lame(10000000, 0.399)
    {
    }

    // Various ways of defining stiffness:
    T bulk_modulus() const { return m_bulk_mod; }
    T poisson() const { return m_poisson; }
    T youngs() const { return m_youngs; }
    T mu() const { return m_mu; }
    T lambda() const { return m_lambda; }
    int model() const { return m_model; }
    T limit_min() const { return m_limit_min; }
    T limit_max() const { return m_limit_max; }

    // Set params with mu and lambda
    void set_lame(T mu_, T lambda_)
    {
        m_mu = mu_;
        m_lambda = lambda_;
        m_youngs = mu_ * (3.0 * lambda_ + 2.0 * mu_) / (mu_ + lambda_);
        m_poisson = lambda_ / (2.0 * (lambda_ + mu_));
        m_bulk_mod = m_lambda + (2.0 / 3.0) * m_mu;
    }

    // Set params with youngs and poisson ratio
    void set_yp(T youngs_, T poisson_)
    {
        m_youngs = youngs_;
        m_poisson = poisson_;
        m_mu = m_youngs / (2.0 * (1.0 + m_poisson));
        m_lambda = m_youngs * m_poisson / ((1.0 + m_poisson) * (1.0 - 2.0 * m_poisson));
        m_bulk_mod = m_lambda + (2.0 / 3.0) * m_mu;
    }

    // Set the elastic model (from the enums above)
    void set_model(int model_) { m_model = model_; }

    // Strain limiting, either hard or soft depending on the model used.
    void set_limits(T min, T max)
    {
        m_limit_min = min;
        m_limit_max = max;
    }

  protected:
    // Hard strain limiting (e.g. [0.95,1.05]), default no limit
    // with  min: -100 to 1, max: 1 to 100.
    // In practice if max>99 it's basically no limiting.
    T m_limit_min, m_limit_max;
    T m_youngs, m_poisson, m_mu, m_lambda, m_bulk_mod;
    int m_model;
};

} // ns mcl
#endif