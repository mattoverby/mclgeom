// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_BARYCOORDS_HPP
#define MCL_BARYCOORDS_HPP 1

#include <Eigen/Core>

namespace mcl {

// vertex-edge barycoords, allows negative barys
template<typename T, int DIM>
static inline Eigen::Matrix<T, 2, 1>
point_edge_barys(const Eigen::Matrix<T, DIM, 1>& p,
                 const Eigen::Matrix<T, DIM, 1>& p0,
                 const Eigen::Matrix<T, DIM, 1>& p1)
{
    Eigen::Matrix<T, DIM, 1> v0 = p1 - p0;
    Eigen::Matrix<T, DIM, 1> v2 = p - p0;
    T d00 = v0.dot(v0);
    if (d00 <= T(0)) {
        return Eigen::Matrix<T, 2, 1>(-1, 0);
    }
    T d20 = v2.dot(v0);
    T invDenom = 1.0 / d00;
    Eigen::Matrix<T, 2, 1> r;
    r[0] = (d00 - d20) * invDenom;
    r[1] = 1.0 - r[0];
    return r;
}

template<typename T, int DIM> // Compute barycentric coords for a point on a triangle
static inline Eigen::Matrix<T, 3, 1>
point_triangle_barys(const Eigen::Matrix<T, DIM, 1>& p,
                     const Eigen::Matrix<T, DIM, 1>& p0,
                     const Eigen::Matrix<T, DIM, 1>& p1,
                     const Eigen::Matrix<T, DIM, 1>& p2)
{
    Eigen::Matrix<T, DIM, 1> v0 = p1 - p0, v1 = p2 - p0, v2 = p - p0;
    T d00 = v0.dot(v0);
    T d01 = v0.dot(v1);
    T d11 = v1.dot(v1);
    T d20 = v2.dot(v0);
    T d21 = v2.dot(v1);
    T invDenom = 1.0 / (d00 * d11 - d01 * d01);
    Eigen::Matrix<T, 3, 1> r;
    r[1] = (d11 * d20 - d01 * d21) * invDenom;
    r[2] = (d00 * d21 - d01 * d20) * invDenom;
    r[0] = 1.0 - r[1] - r[2];
    return r;
}

template<typename T> // Compute barycentric coords for a point in a tet (3D)
static inline Eigen::Matrix<T, 4, 1>
point_tet_barys(const Eigen::Matrix<T, 3, 1>& p,
                const Eigen::Matrix<T, 3, 1>& a,
                const Eigen::Matrix<T, 3, 1>& b,
                const Eigen::Matrix<T, 3, 1>& c,
                const Eigen::Matrix<T, 3, 1>& d)
{
    auto scalar_triple_product = [](const Eigen::Matrix<T, 3, 1>& u,
                                    const Eigen::Matrix<T, 3, 1>& v,
                                    const Eigen::Matrix<T, 3, 1>& w) { return u.dot(v.cross(w)); };
    Eigen::Matrix<T, 3, 1> vap = p - a;
    Eigen::Matrix<T, 3, 1> vbp = p - b;
    Eigen::Matrix<T, 3, 1> vab = b - a;
    Eigen::Matrix<T, 3, 1> vac = c - a;
    Eigen::Matrix<T, 3, 1> vad = d - a;
    Eigen::Matrix<T, 3, 1> vbc = c - b;
    Eigen::Matrix<T, 3, 1> vbd = d - b;
    T va6 = scalar_triple_product(vbp, vbd, vbc);
    T vb6 = scalar_triple_product(vap, vac, vad);
    T vc6 = scalar_triple_product(vap, vad, vab);
    T vd6 = scalar_triple_product(vap, vab, vac);
    T v6 = 1.0 / scalar_triple_product(vab, vac, vad);
    return Eigen::Matrix<T, 4, 1>(va6 * v6, vb6 * v6, vc6 * v6, vd6 * v6);
}

} // end namespace mcl

#endif
