// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_CONSTRAINT_POLISH_HPP
#define MCL_GEOM_CONSTRAINT_POLISH_HPP 1

#include "SignedMeasure.hpp"

#include <Eigen/Core>

namespace mcl {

/// @brief A volume constraint, volume(x) >= target
/// @tparam T scalar type
/// @tparam DIM dimension of vertices, with primitive dimension DIM + 1
template<typename T, int DIM>
class VolumeConstraint
{
public:
    Eigen::Vector<int, DIM+1> stencil = Eigen::Vector<int, DIM+1>::Zero();
    T scaling = 1;

    /// @brief Constructor
    VolumeConstraint() = default;

    /// @brief Constructor, scaling is computed as the surface area if not positive
    VolumeConstraint(const T* x, const Eigen::Vector<int, DIM+1> stencil_, T scaling_ = -1);

    /// @brief Destructor
    ~VolumeConstraint() = default;

    /// @brief Returns nonlinear eval at x = x0*(t-1) + x1*t
    T eval(const T* x0, const T *x1, T t) const;

    /// @brief Returns gradients at x = x0*(t-1) + x1*t
    std::array<Eigen::Vector<T,DIM>, DIM+1> gradients(const T* x0, const T *x1, T t) const;
};

/// @brief For containing and organizing inequality constraints
/// @tparam T scalar type
/// @tparam DIM dimension of vertices, with primitive dimension DIM + 1
template<typename T, int DIM>
class InequalityConstraintSet
{
public:
    std::vector<VolumeConstraint<T,DIM>> volume_constraints;

    /// @brief Clears existing data
    void clear()
    {
        volume_constraints.clear();
    }

    /// @brief Loops over all primitives and adds any that are inverted to the set
    void add_inversions(const T* x, const T* x_rest, const int* primitives, int num_primitives)
    {
        T threshold = T(1e-4);
        volume_constraints.reserve(num_primitives/4);
        for (int i=0; i<num_primitives; ++i) {
            auto stencil = get_primitive<DIM+1>(i, primitives);
            auto verts = get_verts<T,DIM,DIM+1>(x, x, stencil.data(), 0);
            if constexpr (DIM == 2) {
                if (signed_triangle_area(verts[0], verts[1], verts[2]) <= threshold) {
                    volume_constraints.emplace_back(x_rest, stencil);
                }
            }
            else if constexpr (DIM == 3) {
                if (signed_tet_volume(verts[0], verts[1], verts[2], verts[3]) <= threshold) {
                    volume_constraints.emplace_back(x_rest, stencil);
                }
            }
        }
    }
};

//
// Implementation
//

template<typename T, int DIM>
VolumeConstraint<T,DIM>::VolumeConstraint(const T* x, const Eigen::Vector<int, DIM+1> stencil_, T scaling_) : stencil(stencil_), scaling(scaling_)
{
    if (scaling <= T(0)) {
        auto verts = get_verts<T,DIM,DIM+1>(x, x, stencil.data(), 0);
        if constexpr (DIM == 2) {
            scaling = T(1) / (T(0.5) * triangle_perimeter(verts[0], verts[1], verts[2]));
        }
        else if constexpr (DIM == 3) {
            scaling =  T(1) / (T(0.5) * tet_surface_area(verts[0], verts[1], verts[2], verts[2]));
        }
    }
}

template<typename T, int DIM>
T VolumeConstraint<T,DIM>::eval(const T* x0, const T *x1, T t) const
{
    if constexpr (DIM == 2)
    {
        auto verts = get_verts<T,DIM,3>(0, x0, x1, stencil.data(), t);
        return scaling * signed_triangle_area(verts[0], verts[1], verts[2]);
    }
    else if constexpr (DIM == 3)
    {
        auto verts = get_verts<T,DIM,4>(0, x0, x1, stencil.data(), t);
        return scaling * signed_tet_volume(verts[0], verts[1], verts[2], verts[2]);
    }
    return 0;  
}

template<typename T, int DIM>
std::array<Eigen::Vector<T,DIM>, DIM+1> VolumeConstraint<T,DIM>::gradients(const T* x0, const T *x1, T t) const
{
    if constexpr (DIM == 2)
    {
        auto verts = get_verts<T,DIM,3>(0, x0, x1, stencil.data(), t);
        auto grads = signed_triangle_area_gradients(verts[0], verts[1], verts[2]);
        grads[0] *= scaling;
        grads[1] *= scaling;
        grads[2] *= scaling;
        return grads;
    }
    else if constexpr (DIM == 3)
    {
        auto verts = get_verts<T,DIM,4>(0, x0, x1, stencil.data(), t);
        auto grads = signed_tet_volume_gradients(verts[0], verts[1], verts[2], verts[2]);
        grads[0] *= scaling;
        grads[1] *= scaling;
        grads[2] *= scaling;
        grads[3] *= scaling;
        return grads;
    }
    return {{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}};
}


} // end ns mcl

#endif // MCL_GEOM_CONSTRAINT_POLISH_HPP