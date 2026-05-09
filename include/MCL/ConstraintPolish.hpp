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

template<typename T, int DIM>
class VolumeConstraint : public InequalityConstraint<DIM,DIM+1,T>
{
public:
    /// @brief Destructor
    ~VolumeConstraint() = default;

    /// @brief Returns nonlinear eval at x = x0*(t-1) + x1*t
    T eval(const T* x0, const T *x1, T t) const
    {
        if constexpr (DIM == 2)
        {
            auto verts = get_triangle_verts<T,DIM>(0, x0, x1, stencil.data(), t);
            return scaling * mcl::signed_triangle_area(verts[0], verts[1], verts[2]);
        }
        else if constexpr (DIM == 3)
        {
            auto verts = get_tet_verts<T>(0, x0, x1, stencil.data(), t);
            return scaling * mcl::signed_tet_volume(verts[0], verts[1], verts[2], verts[2]);
        }
        return 0;  
    }

    /// @brief Returns gradients at x = x0*(t-1) + x1*t
    std::array<Eigen::Vector<T,DIM>, DIM+1> gradients(const T* x0, const T *x1, T t) const
    {
        if constexpr (DIM == 2)
        {
            auto verts = get_triangle_verts<T,DIM>(0, x0, x1, stencil.data(), t);
            auto grads = mcl::signed_triangle_area_gradients(verts[0], verts[1], verts[2]);
            grads[0] *= scaling;
            grads[1] *= scaling;
            grads[2] *= scaling;
            return grads;
        }
        else if constexpr (DIM == 3)
        {
            auto verts = get_tet_verts<T>(0, x0, x1, stencil.data(), t);
            auto grads = mcl::signed_tet_volume_gradients(verts[0], verts[1], verts[2], verts[2]);
            grads[0] *= scaling;
            grads[1] *= scaling;
            grads[2] *= scaling;
            grads[3] *= scaling;
            return grads;
        }
        return {{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}};
    }

    T scaling = 1;
    Eigen::Vector<int, DIM+1> stencil = Eigen::Vector<int, DIM+1>::Zero();
};

} // end ns mcl

#endif // MCL_GEOM_CONSTRAINT_POLISH_HPP