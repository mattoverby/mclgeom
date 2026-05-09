// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_CONSTRAINT_POLISH_HPP
#define MCL_GEOM_CONSTRAINT_POLISH_HPP 1

#include "SignedMeasure.hpp"

#include <Eigen/Core>

namespace mcl {

template<int DIM, int SIZE, typename T>
class InequalityConstraint
{
public:
    /// @brief Destructor
    virtual ~InequalityConstraint() = default;

    /// @brief Returns nonlinear eval at x = x0*(t-1) + x1*t
    virtual T eval(const T* x0, const T *x1, T t) const = 0;

    /// @brief Returns gradients at x = x0*(t-1) + x1*t
    std::array<Eigen::Vector<T,DIM>, SIZE> gradients(const T* x0, const T *x1, T t) const = 0;
};

template<typename T>
class VolumeConstraint : public InequalityConstraint<3,4,T>
{
public:
    /// @brief Destructor
    ~VolumeConstraint() = default;

    /// @brief Returns nonlinear eval at x = x0*(t-1) + x1*t
    T eval(const T* x0, const T *x1, T t) const
    {
        auto verts = get_tet_verts(0, x0, x1, stencil.data(), 0);

    }

    /// @brief Returns gradients at x = x0*(t-1) + x1*t
    std::array<Eigen::Vector3<T>, 4> gradients(const T* x0, const T *x1, T t) const
    {
        auto verts = get_tet_verts(0, x0, x1, stencil.data(), 0);

    }


    Eigen::Vector4i stencil;
};

} // end ns mcl

#endif // MCL_GEOM_CONSTRAINT_POLISH_HPP