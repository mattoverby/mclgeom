// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_NORMAL_HPP
#define MCL_GEOM_NORMAL_HPP 1

#include <Eigen/Dense>

namespace mcl {

/// @brief Returns the triangle normal
template<typename T>
inline Eigen::Vector3<T>
triangle_normal(const Eigen::Vector3<T>& a,
                const Eigen::Vector3<T>& b,
                const Eigen::Vector3<T>& c,
                bool unit_length = true)
{
    Eigen::Vector3<T> n = (b - a).cross(c - a);
    if (unit_length) {
        T len = n.norm();
        if (len > 0) {
            n /= len;
        }
    }
    return n;
}

/// @brief Returns the 2D edge normal
template<typename T>
inline Eigen::Vector2<T>
edge_normal(const Eigen::Vector2<T>& p0, const Eigen::Vector2<T>& p1, bool unit_length = true)
{
    Eigen::Vector2<T> n(p1[1] - p0[1], -(p1[0] - p0[0]));
    if (unit_length) {
        T len = n.norm();
        if (len > 0) {
            n /= len;
        }
    }
    return n;
}

} // ns mcl

#endif
