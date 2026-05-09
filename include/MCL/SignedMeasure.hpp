// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_SIGNEDMEASURE_HPP
#define MCL_SIGNEDMEASURE_HPP 1

#include <Eigen/Dense>

namespace mcl {

/// @brief Returns the signed area of a 2D triangle
template <typename T>
inline T
signed_triangle_area(const Eigen::Vector2<T>& p1, const Eigen::Vector2<T>& p2, const Eigen::Vector2<T>& p3);

/// @brief Returns (unscaled) first derivative of 2D signed triangle area
template <typename T>
inline std::array<Eigen::Vector2<T>, 3>
signed_triangle_area_gradients(const Eigen::Vector2<T>& p1, const Eigen::Vector2<T>& p2, const Eigen::Vector2<T>& p3);

inline double
triangle_perimeter(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3);

inline double
triangle_area(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3);

inline double
signed_tet_volume(const Eigen::Vector3d& p1,
                  const Eigen::Vector3d& p2,
                  const Eigen::Vector3d& p3,
                  const Eigen::Vector3d& p4);

// Returns (unscaled) first derivative of signed tet volume
inline std::array<Eigen::Vector3d, 4>
signed_tet_volume_gradients(const Eigen::Vector3d& p1,
                            const Eigen::Vector3d& p2,
                            const Eigen::Vector3d& p3,
                            const Eigen::Vector3d& p4);
inline double
tet_surface_area(const Eigen::Vector3d& p1,
                 const Eigen::Vector3d& p2,
                 const Eigen::Vector3d& p3,
                 const Eigen::Vector3d& p4);

/// @brief Helper accessor returns tet vertices x = x0*(1-t) + x1*t
template <typename T, typename IndexType>
inline std::array<Eigen::Vector3<T>, 4>
get_tet_verts(int tet_index, const T *V0, const T *V1, const IndexType *tets, T t);

//
// Implementation
//

template <typename T>
T
signed_triangle_area(const Eigen::Vector2<T>& p1, const Eigen::Vector2<T>& p2, const Eigen::Vector2<T>& p3)
{
    return T(0.5) * (-p2[0] * p1[1] + p3[0] * p1[1] + p1[0] * p2[1] - p3[0] * p2[1] - p1[0] * p3[1] + p2[0] * p3[1]);
}

template <typename T>
std::array<Eigen::Vector2<T>, 3>
signed_triangle_area_gradients(const Eigen::Vector2<T>& a, const Eigen::Vector2<T>& b, const Eigen::Vector2<T>& c)
{
    std::array<Eigen::Vector2<T>, 3> g;
    g[0] = T(0.5) * Eigen::Vector2<T>(b[1] - c[1], -b[0] + c[0]);
    g[1] = T(0.5) * Eigen::Vector2<T>(-a[1] + c[1], a[0] - c[0]);
    g[2] = T(0.5) * Eigen::Vector2<T>(a[1] - b[1], -a[0] + b[0]);
    return g;
}

inline double
triangle_perimeter(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const Eigen::Vector2d& p3)
{
    return (p1 - p2).norm() + (p2 - p3).norm() + (p3 - p1).norm();
}

// https://en.wikipedia.org/wiki/Heron%27s_formula
inline double
triangle_area(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3)
{
    double a = (p1 - p2).norm();
    double b = (p2 - p3).norm();
    double c = (p3 - p1).norm();
    double s = (a + b + c) * 0.5;
    return std::sqrt(s * (s - a) * (s - b) * (s - c));
}

inline double
signed_tet_volume(const Eigen::Vector3d& p1,
                  const Eigen::Vector3d& p2,
                  const Eigen::Vector3d& p3,
                  const Eigen::Vector3d& p4)
{
    Eigen::Matrix3d edges;
    edges.col(0) = p2 - p1;
    edges.col(1) = p3 - p1;
    edges.col(2) = p4 - p1;
    return (1.0 / 6.0) * edges.determinant();
}

inline std::array<Eigen::Vector3d, 4>
signed_tet_volume_gradients(const Eigen::Vector3d& a,
                            const Eigen::Vector3d& b,
                            const Eigen::Vector3d& c,
                            const Eigen::Vector3d& d)
{
    std::array<Eigen::Vector3d, 4> grads;
    const Eigen::Vector3d& p0 = a;
    const Eigen::Vector3d& p1 = b;
    const Eigen::Vector3d& p2 = c;
    const Eigen::Vector3d& p3 = d;
    static const double sixth = (1.0 / 6.0);
    grads[0] = sixth * (p1 - p2).cross(p3 - p2);
    grads[1] = sixth * (p2 - p0).cross(p3 - p0);
    grads[2] = sixth * (p0 - p1).cross(p3 - p1);
    grads[3] = sixth * (p1 - p0).cross(p2 - p0);
    return grads;
}

inline double
tet_surface_area(const Eigen::Vector3d& p1,
                 const Eigen::Vector3d& p2,
                 const Eigen::Vector3d& p3,
                 const Eigen::Vector3d& p4)
{
    double a1 = triangle_area(p2, p3, p4);
    double a2 = triangle_area(p2, p3, p1);
    double a3 = triangle_area(p3, p4, p1);
    double a4 = triangle_area(p4, p2, p1);
    return (a1 + a2 + a3 + a4);
}

template <typename T, typename IndexType>
std::array<Eigen::Vector3<T>, 4>
get_tet_verts(int tet_index, const T *V0, const T *V1, const IndexType *tets, T t)
{
    Eigen::Vector4i tet(tets[tet_index*4+0], tets[tet_index*4+1], tets[tet_index*4+2], tets[tet_index*4+3]);
    std::array<Eigen::Vector3<T>, 4> v0 = {
        Eigen::Vector3<T>(V0[tet[0]*3+0], V0[tet[0]*3+1], V0[tet[0]*3+2]),
        Eigen::Vector3<T>(V0[tet[1]*3+0], V0[tet[1]*3+1], V0[tet[1]*3+2]),
        Eigen::Vector3<T>(V0[tet[2]*3+0], V0[tet[2]*3+1], V0[tet[2]*3+2]),
        Eigen::Vector3<T>(V0[tet[3]*3+0], V0[tet[3]*3+1], V0[tet[3]*3+2])};
    std::array<Eigen::Vector3<T>, 4> v1 = {
        Eigen::Vector3<T>(V1[tet[0]*3+0], V1[tet[0]*3+1], V1[tet[0]*3+2]),
        Eigen::Vector3<T>(V1[tet[1]*3+0], V1[tet[1]*3+1], V1[tet[1]*3+2]),
        Eigen::Vector3<T>(V1[tet[2]*3+0], V1[tet[2]*3+1], V1[tet[2]*3+2]),
        Eigen::Vector3<T>(V1[tet[3]*3+0], V1[tet[3]*3+1], V1[tet[3]*3+2])};
        return {
            v0[0]*(T(1) - t) + v1[0]*t,
            v0[1]*(T(1) - t) + v1[1]*t,
            v0[2]*(T(1) - t) + v1[2]*t,
            v0[3]*(T(1) - t) + v1[3]*t };
}

} // ns mcl

#endif
