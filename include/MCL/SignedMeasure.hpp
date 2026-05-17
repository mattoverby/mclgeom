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

/// @brief Returns the 2D triangle perimeter
template <typename T>
inline T
triangle_perimeter(const Eigen::Vector2<T>& p1, const Eigen::Vector2<T>& p2, const Eigen::Vector2<T>& p3);

/// @brief Returns the 2D triangle area
template <typename T>
inline T
triangle_area(const Eigen::Vector3<T>& p1, const Eigen::Vector3<T>& p2, const Eigen::Vector3<T>& p3);

/// @brief Returns signed volume of the tet
template <typename T>
inline T
signed_tet_volume(const Eigen::Vector3<T>& p1,
                  const Eigen::Vector3<T>& p2,
                  const Eigen::Vector3<T>& p3,
                  const Eigen::Vector3<T>& p4);

/// @brief Returns (unscaled) first derivative of signed tet volume
template <typename T>
inline std::array<Eigen::Vector3<T>, 4>
signed_tet_volume_gradients(const Eigen::Vector3<T>& p1,
                            const Eigen::Vector3<T>& p2,
                            const Eigen::Vector3<T>& p3,
                            const Eigen::Vector3<T>& p4);

/// @brief Returns the tet surface area                            
template <typename T>
inline T
tet_surface_area(const Eigen::Vector3<T>& p1,
                 const Eigen::Vector3<T>& p2,
                 const Eigen::Vector3<T>& p3,
                 const Eigen::Vector3<T>& p4);

/// @brief Helper accessor returns vertex refs
template <typename T, int DIM, int PDIM>
inline std::array<Eigen::Vector<T, DIM>, PDIM>
get_verts(const T *V, const int *stencil);

/// @brief Helper accessor returns vertices x = x0*(1-t) + x1*t
template <typename T, int DIM, int PDIM>
inline std::array<Eigen::Vector<T,DIM>, PDIM>
get_verts_at_delta(const T *V0, const T *V1, const int *stencil, T t);

/// @brief Helper accessor returns the primitive
template <int PDIM>
inline Eigen::Vector<int, PDIM>
get_primitive(int prim_index, const int *primitives);


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

template <typename T>
T
triangle_perimeter(const Eigen::Vector2<T>& p1, const Eigen::Vector2<T>& p2, const Eigen::Vector2<T>& p3)
{
    return (p1 - p2).norm() + (p2 - p3).norm() + (p3 - p1).norm();
}

// https://en.wikipedia.org/wiki/Heron%27s_formula
template <typename T>
T
triangle_area(const Eigen::Vector3<T>& p1, const Eigen::Vector3<T>& p2, const Eigen::Vector3<T>& p3)
{
    T a = (p1 - p2).norm();
    T b = (p2 - p3).norm();
    T c = (p3 - p1).norm();
    T s = (a + b + c) * 0.5;
    return std::sqrt(s * (s - a) * (s - b) * (s - c));
}

template <typename T>
T
signed_tet_volume(const Eigen::Vector3<T>& p0,
                  const Eigen::Vector3<T>& p1,
                  const Eigen::Vector3<T>& p2,
                  const Eigen::Vector3<T>& p3)
{
    Eigen::Matrix3<T> edges;
    edges.col(0) = p1 - p0;
    edges.col(1) = p2 - p0;
    edges.col(2) = p3 - p0;
    return T(1.0 / 6.0) * edges.determinant();
}

template <typename T>
std::array<Eigen::Vector3<T>, 4>
signed_tet_volume_gradients(const Eigen::Vector3<T>& p0,
                            const Eigen::Vector3<T>& p1,
                            const Eigen::Vector3<T>& p2,
                            const Eigen::Vector3<T>& p3)
{
    std::array<Eigen::Vector3<T>, 4> grads;
    static const T sixth = T(1.0 / 6.0);
    grads[0] = sixth * (p1 - p2).cross(p3 - p2);
    grads[1] = sixth * (p2 - p0).cross(p3 - p0);
    grads[2] = sixth * (p0 - p1).cross(p3 - p1);
    grads[3] = sixth * (p1 - p0).cross(p2 - p0);
    return grads;
}

template <typename T>
T
tet_surface_area(const Eigen::Vector3<T>& p1,
                 const Eigen::Vector3<T>& p2,
                 const Eigen::Vector3<T>& p3,
                 const Eigen::Vector3<T>& p4)
{
    T a1 = triangle_area(p2, p3, p4);
    T a2 = triangle_area(p2, p3, p1);
    T a3 = triangle_area(p3, p4, p1);
    T a4 = triangle_area(p4, p2, p1);
    return (a1 + a2 + a3 + a4);
}

template <typename T, int DIM, int PDIM>
std::array<Eigen::Vector<T, DIM>, PDIM>
get_verts(const T *V, const int *stencil)
{
    std::array<Eigen::Vector<T, DIM>, PDIM> v;
    for (int i = 0; i < PDIM; ++i) {
        v[i] = Eigen::Map<const Eigen::Vector<T, DIM>>(V + stencil[i] * DIM).eval();
    }
    return v;
}

template <typename T, int DIM, int PDIM>
std::array<Eigen::Vector<T,DIM>, PDIM>
get_verts_at_delta(const T *V0, const T *V1, const int *stencil, T t)
{
    auto verts_t0 = get_verts<T,DIM,PDIM>(V0, stencil);
    auto verts_t1 = get_verts<T,DIM,PDIM>(V1, stencil);
    std::array<Eigen::Vector<T,DIM>, PDIM> v;
    for (int i = 0; i < PDIM; ++i) {
        v[i] = verts_t0[i] * (T(1) - t) + verts_t1[i] * t;
    }
    return v;
}

template <int PDIM>
Eigen::Vector<int, PDIM>
get_primitive(int prim_index, const int *primitives)
{
    return Eigen::Map<const Eigen::Vector<int, PDIM>>(primitives + prim_index * PDIM).eval();
}

} // ns mcl

#endif
