// Copyright Matt Overby 2025.
// Distributed under the MIT License.

#ifndef MCL_GEOM_BENDINGMODEL_HPP
#define MCL_GEOM_BENDINGMODEL_HPP 1

#include <Eigen/Dense>
#include <map>

namespace mcl {

/// @brief Given triangle F, computes hinges H, with each row of H a four-stencil:
/// (x0, x1) is the shared edge between triangles (x0, x1, x2) and (x1, x0, x3).
/// Triangles on the boundary are ommitted, or optionally have x2 == x3.
/// Assumes but does not check for non-manifold mesh.
template<typename DerivedF, typename DerivedH>
inline void
make_hinges(const DerivedF& F, DerivedH& H, bool include_boundary = false);

/// @brief Computes matrix "Q" of DOI 10.5555/1281957.1281987
/// "A quadratic bending model for inextensible surfaces", Bergou et al.
/// (x0, x1) is the shared edge between triangles (x0, x1, x2) and (x1, x0, x3).
template<typename T>
inline Eigen::Matrix4<T>
quadratic_bend_Q(const Eigen::Vector3<T>& x0,
                 const Eigen::Vector3<T>& x1,
                 const Eigen::Vector3<T>& x2,
                 const Eigen::Vector3<T>& x3);

/// @brief Computes row-vector "K0" of DOI 10.5555/1281957.1281987
/// "A quadratic bending model for inextensible surfaces", Bergou et al.
/// (x0, x1) is the shared edge between triangles (x0, x1, x2) and (x1, x0, x3).
template<typename T>
inline Eigen::RowVector4<T>
quadratic_bend_K0(const Eigen::Vector3<T>& x0,
                  const Eigen::Vector3<T>& x1,
                  const Eigen::Vector3<T>& x2,
                  const Eigen::Vector3<T>& x3);

/// @brief Computes alpha coeffs of DOI 10.5555/1218064.1218078
/// "Simple Linear Bending Stiffness in Particle Systems", Volino & Magnenat-Thalmann.
/// (Pc, Pd) is the shared edge between triangles (Pa, Pd, Pc) and (Pb, Pc, Pd).
/// Note: this is different than quadratic model and "make_hinges" above!
template<typename T>
Eigen::Vector4<T>
linear_bend_alpha(const Eigen::Vector3<T>& Pa,
                  const Eigen::Vector3<T>& Pb,
                  const Eigen::Vector3<T>& Pc,
                  const Eigen::Vector3<T>& Pd);

//
// Implementation
//

template<typename DerivedF, typename DerivedH>
void
make_hinges(const DerivedF& F, DerivedH& H, bool include_boundary)
{
    struct Entry
    {
        int x0, x1;
        int opp;
    };

    std::map<std::pair<int, int>, Entry> edge_map;
    std::vector<Eigen::RowVector4i> hinges;
    hinges.reserve(F.rows());

    for (int f = 0; f < int(F.rows()); ++f) {
        int v[3] = { F(f, 0), F(f, 1), F(f, 2) };

        for (int e = 0; e < 3; ++e) {
            int x0 = v[e];
            int x1 = v[(e + 1) % 3];
            int x2 = v[(e + 2) % 3];
            std::pair<int, int> key = std::make_pair(x0, x1);
            if (x1 > x0) {
                std::swap(key.first, key.second);
            }

            // See if the edge exists. If not, add it.
            auto it = edge_map.find(key);
            if (it == edge_map.end()) {
                edge_map.emplace(key, Entry{ x0, x1, x2 });
            }
            // If it exists, form the hinge
            else {
                const Entry& e0 = it->second;
                hinges.emplace_back(e0.x0, e0.x1, e0.opp, x2);
                edge_map.erase(it);
            }
        }
    }

    // boundary edges, what hasn't been removed in above loop
    if (include_boundary) {
        for (const auto& kv : edge_map) {
            const Entry& e = kv.second;
            hinges.emplace_back(e.x0, e.x1, e.opp, e.opp);
        }
    }

    H.resize(hinges.size(), 4);
    for (int i = 0; i < int(hinges.size()); ++i) {
        H.row(i) = hinges[i];
    }
}

template<typename T>
Eigen::Matrix4<T>
quadratic_bend_Q(const Eigen::Vector3<T>& x0,
                 const Eigen::Vector3<T>& x1,
                 const Eigen::Vector3<T>& x2,
                 const Eigen::Vector3<T>& x3)
{
    Eigen::RowVector4<T> K0 = quadratic_bend_K0(x0, x1, x2, x3);
    T A1 = 0.5 * ((x1 - x0).cross(x2 - x0)).norm();
    T A2 = 0.5 * ((x0 - x1).cross(x3 - x1)).norm();
    T scale = std::abs(A1 + A2) > 0 ? T(3 * (A1 + A2)) : T(1);
    Eigen::Matrix4<T> Q = scale * (K0.transpose() * K0);
    return Q;
}

template<typename T>
Eigen::RowVector4<T>
quadratic_bend_K0(const Eigen::Vector3<T>& x0,
                  const Eigen::Vector3<T>& x1,
                  const Eigen::Vector3<T>& x2,
                  const Eigen::Vector3<T>& x3)
{
    // Computes cotangent between (v1-v0) and (v2-v0)
    auto cotangent = [&](const Eigen::Vector3<T>& v0, const Eigen::Vector3<T>& v1, const Eigen::Vector3<T>& v2) {
        Eigen::Vector3<T> e0 = v1 - v0;
        Eigen::Vector3<T> e1 = v2 - v0;
        T cos_theta = e0.dot(e1);
        T sin_theta = e0.cross(e1).norm();
        return std::abs(sin_theta) > 0 ? cos_theta / sin_theta : T(0);
    };

    T c01 = cotangent(x0, x1, x2);
    T c03 = cotangent(x1, x0, x2);
    T c02 = cotangent(x0, x3, x1);
    T c04 = cotangent(x1, x3, x0);
    Eigen::RowVector4<T> K0;
    K0[0] = (c03 + c04);
    K0[1] = (c01 + c02);
    K0[2] = -(c01 + c03);
    K0[3] = -(c02 + c04);
    return K0;
}

template<typename T>
Eigen::Vector4<T>
linear_bend_alpha(const Eigen::Vector3<T>& Pa,
                  const Eigen::Vector3<T>& Pb,
                  const Eigen::Vector3<T>& Pc,
                  const Eigen::Vector3<T>& Pd)
{
    // Normals: Eq. 6
    Eigen::Vector3<T> NA = (Pa - Pc).cross(Pa - Pd);
    Eigen::Vector3<T> NB = (Pb - Pd).cross(Pb - Pc);
    Eigen::Vector3<T> NC = (Pc - Pb).cross(Pc - Pa);
    Eigen::Vector3<T> ND = (Pd - Pa).cross(Pd - Pb);
    T lenNA = NA.norm();
    T lenNB = NB.norm();
    T lenNC = NC.norm();
    T lenND = ND.norm();

    if (lenNA == T(0) || lenNB == T(0) || lenNC == T(0) || lenND == T(0)) {
        return Eigen::Vector4d::Zero();
    }

    // Coeffs: Eq. 5
    T hA = lenNA;
    T hB = lenNB;
    T sumN = lenNA + lenNB;
    T sumC = lenNC + lenND;
    T alphaA = lenNB / (sumN + sumC);
    T alphaB = -lenNA / (sumN + sumC);
    T alphaC = lenNB / (sumN + sumC);
    T alphaD = -lenNA / (sumN + sumC);

    return Eigen::Vector4d(alphaA, alphaB, alphaC, alphaD);
}

} // end namespace mcl

#endif // MCL_GEOM_BENDINGMODEL_HPP