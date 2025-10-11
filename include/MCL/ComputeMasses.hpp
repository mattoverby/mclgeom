// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_COMPUTEMASSES_HPP
#define MCL_COMPUTEMASSES_HPP 1

#include <Eigen/Dense>

namespace mcl {

// V are vertices (n x 2 or 3)
// P are primitives (m x 3 or 4)
// M is n x 1 of per-vertex masses
// Computes volume (or area) weighted masses with unit-volume density.
// If negative, defaults are used: 1100 for volumetric, 0.4 for cloth/2D
// See: https://www.engineeringtoolbox.com/density-solids-d_1265.html.
// Masses between the min and max indices in P are set.
// If unreferenced, they are set to zero.
// Returns true if all masses in span are positive.
template<typename DerivedV, typename DerivedP, typename Scalar>
static inline bool
compute_masses(const Eigen::MatrixBase<DerivedV>& V,
               const Eigen::MatrixBase<DerivedP>& P,
               Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& M,
               double density_kgd = -1)
{
    using namespace Eigen;
    int V_dim = V.cols();
    int P_dim = P.cols();
    if (V_dim < 2 || V_dim > 3) {
        return false;
    }
    if (P_dim < 3 || P_dim > 4) {
        return false;
    }

    // Use 3D vec for calculation even if 2D
    auto Vi = [&](int idx) {
        Vector3d vi = Vector3d::Zero();
        vi.head(V_dim) = V.row(idx).template cast<double>();
        return vi;
    };

    // Default densities
    if (density_kgd < 0) {
        if (V_dim == 2 || P_dim == 3) {
            density_kgd = 0.4;
        } else if (V_dim == 3 && P_dim == 4) {
            density_kgd = 1100;
        }
    }

    // Resize masses if needed and set span to zero
    int min_Pi = P.minCoeff();
    int max_Pi = P.maxCoeff();
    if (M.rows() < max_Pi) {
        M.conservativeResize(max_Pi + 1);
    }
    M.segment(min_Pi, max_Pi - min_Pi + 1).array() = 0;

    // Compute mass contrib from each element
    int np = P.rows();
    for (int i = 0; i < np; ++i) {
        if (P_dim == 4) {
            Vector3d p_verts[4] = { Vi(P(i, 0)), Vi(P(i, 1)), Vi(P(i, 2)), Vi(P(i, 3)) };
            Matrix<double, 3, 3> E;
            E.col(0) = p_verts[1] - p_verts[0];
            E.col(1) = p_verts[2] - p_verts[0];
            E.col(2) = p_verts[3] - p_verts[0];
            double vol = std::abs(E.determinant() / 6.0);
            double tet_mass = density_kgd * vol;
            M[P(i, 0)] += Scalar(tet_mass / 4.0);
            M[P(i, 1)] += Scalar(tet_mass / 4.0);
            M[P(i, 2)] += Scalar(tet_mass / 4.0);
            M[P(i, 3)] += Scalar(tet_mass / 4.0);
        } else if (P_dim == 3) {
            Vector3d p_verts[3] = { Vi(P(i, 0)), Vi(P(i, 1)), Vi(P(i, 2)) };
            Vector3d e0 = p_verts[1] - p_verts[0];
            Vector3d e1 = p_verts[2] - p_verts[0];
            double area = 0.5 * (e0.cross(e1)).norm();
            double tri_mass = density_kgd * area;
            M[P(i, 0)] += Scalar(tri_mass / 3.0);
            M[P(i, 1)] += Scalar(tri_mass / 3.0);
            M[P(i, 2)] += Scalar(tri_mass / 3.0);
        }
    }

    double min_mass = M.segment(min_Pi, max_Pi - min_Pi + 1).minCoeff();
    return min_mass > 0;
}

} // ns mcl

#endif
