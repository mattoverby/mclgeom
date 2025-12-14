// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_COMPUTEMASSES_HPP
#define MCL_GEOM_COMPUTEMASSES_HPP

#include <Eigen/Dense>
#include <limits>

namespace mcl {

// V are vertices (n x 2 or 3)
// P are primitives (m x 2, 3 or 4 for edges/triangles/tets)
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
               Scalar density_kgd = -1)
{
    int V_dim = V.cols();
    int P_dim = P.cols();
    if (V_dim < 2 || V_dim > 3) {
        return false;
    }
    if (P_dim < 2 || P_dim > 4) {
        return false;
    }

    // Use 3D vec for calculation even if 2D
    auto Vi = [&](int idx) {
        Eigen::Vector3<Scalar> vi = Eigen::Vector3<Scalar>::Zero();
        vi.head(V_dim) = V.row(idx).template cast<Scalar>();
        return vi;
    };

    // Default densities
    if (density_kgd < 0) {
        if (P_dim == 2) { // approximate by area/volume
            density_kgd = V_dim == 2 ? 0.4 : 1000.0;
        } else if (V_dim == 2 || P_dim == 3) {
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

    // If P_dim == 2, approximate by area/volume
    Eigen::Vector3<Scalar> bmin, bmax;
    for (int i = 0; i < 3; ++i) {
        bmin[i] = std::numeric_limits<Scalar>::max();
        bmax[i] = std::numeric_limits<Scalar>::lowest();
    }

    // Compute mass contrib from each element
    int np = P.rows();
    for (int i = 0; i < np; ++i) {
        if (P_dim == 4) {
            Eigen::Vector3<Scalar> p_verts[4] = { Vi(P(i, 0)), Vi(P(i, 1)), Vi(P(i, 2)), Vi(P(i, 3)) };
            Eigen::Matrix<Scalar, 3, 3> E;
            E.col(0) = p_verts[1] - p_verts[0];
            E.col(1) = p_verts[2] - p_verts[0];
            E.col(2) = p_verts[3] - p_verts[0];
            Scalar vol = std::abs(E.determinant() / 6.0);
            Scalar tet_mass = density_kgd * vol;
            M[P(i, 0)] += Scalar(tet_mass / 4.0);
            M[P(i, 1)] += Scalar(tet_mass / 4.0);
            M[P(i, 2)] += Scalar(tet_mass / 4.0);
            M[P(i, 3)] += Scalar(tet_mass / 4.0);
        } else if (P_dim == 3) {
            Eigen::Vector3<Scalar> p_verts[3] = { Vi(P(i, 0)), Vi(P(i, 1)), Vi(P(i, 2)) };
            Eigen::Vector3<Scalar> e0 = p_verts[1] - p_verts[0];
            Eigen::Vector3<Scalar> e1 = p_verts[2] - p_verts[0];
            Scalar area = 0.5 * (e0.cross(e1)).norm();
            Scalar tri_mass = density_kgd * area;
            M[P(i, 0)] += Scalar(tri_mass / 3.0);
            M[P(i, 1)] += Scalar(tri_mass / 3.0);
            M[P(i, 2)] += Scalar(tri_mass / 3.0);
        } else if (P_dim == 2) {
            Eigen::Vector3<Scalar> p_verts[2] = { Vi(P(i, 0)), Vi(P(i, 1)) };
            for (int j = 0; j < 3; ++j) {
                bmin[j] = std::min(bmin[j], std::min(p_verts[0][j], p_verts[1][j]));
                bmax[j] = std::max(bmax[j], std::max(p_verts[0][j], p_verts[1][j]));
            }
        }
    }

    // TODO: Better mass computation for springs
    if (P_dim == 2 && V_dim == 2) {
        Scalar area = (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]);
        M.array() = density_kgd * area / Scalar(M.rows());
    } else if (P_dim == 2 && V_dim == 3) {
        Scalar volume = (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[2] - bmin[2]);
        M.array() = density_kgd * volume / Scalar(M.rows());
    }

    Scalar min_mass = M.segment(min_Pi, max_Pi - min_Pi + 1).minCoeff();
    return min_mass > 0;
}

} // ns mcl

#endif
