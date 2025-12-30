// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_COMPUTEMASSES_HPP
#define MCL_GEOM_COMPUTEMASSES_HPP

#include <Eigen/Dense>
#include <limits>

namespace mcl {

/// @brief Computes area/volume weighted vertex masses w.r.t. density.
/// For density, see: https://www.engineeringtoolbox.com/density-solids-d_1265.html.
/// If a vertex is unreferenced, its mass will be zero.
/// @return true if all masses > 0.
template<typename DerivedV, typename DerivedP, typename Scalar>
inline bool
compute_masses(const Eigen::MatrixBase<DerivedV>& V,
               const Eigen::MatrixBase<DerivedP>& P,
               Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& M,
               Scalar density_kgd = -1)
{
    enum class MassType
    {
        UNKNOWN,
        SPRING, // TODO
        TRIANGLE,
        TET
    };

    // See what domain we are working with
    int V_dim = V.cols();
    int P_dim = P.cols();
    auto mass_type = MassType::UNKNOWN;
    if (P_dim == 2){
        mass_type = MassType::SPRING;
    }
    else if(P_dim == 3) {
         mass_type = MassType::TRIANGLE;
    }
    else if (P_dim == 4){
      mass_type = MassType::TET;
    }

    if (mass_type == MassType::SPRING || mass_type == MassType::UNKNOWN){
        return false;
    }

    // Use 3D vec for calculation even if 2D to simplify
    auto vertex_as_3D = [&](int idx) {
        Eigen::Vector3<Scalar> vi = Eigen::Vector3<Scalar>::Zero();
        vi.head(V_dim) = V.row(idx).template cast<Scalar>();
        return vi;
    };

    // Default densities
    if (density_kgd < 0) {
        if (mass_type == MassType::TRIANGLE) {
            density_kgd = 0.4;
        } else if (mass_type == MassType::TET) {
            density_kgd = 1100;
        }
    }

    M = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero(V.rows());

    // Compute mass contrib from each element
    int np = P.rows();
    for (int i = 0; i < np; ++i) {
        if (mass_type == MassType::TET) {
            Eigen::Vector3<Scalar> p_verts[4] = { vertex_as_3D(P(i, 0)), vertex_as_3D(P(i, 1)), vertex_as_3D(P(i, 2)), vertex_as_3D(P(i, 3)) };
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
        } else if (mass_type == MassType::TRIANGLE) {
            Eigen::Vector3<Scalar> p_verts[3] = { vertex_as_3D(P(i, 0)), vertex_as_3D(P(i, 1)), vertex_as_3D(P(i, 2)) };
            Eigen::Vector3<Scalar> e0 = p_verts[1] - p_verts[0];
            Eigen::Vector3<Scalar> e1 = p_verts[2] - p_verts[0];
            Scalar area = 0.5 * (e0.cross(e1)).norm();
            Scalar tri_mass = density_kgd * area;
            M[P(i, 0)] += Scalar(tri_mass / 3.0);
            M[P(i, 1)] += Scalar(tri_mass / 3.0);
            M[P(i, 2)] += Scalar(tri_mass / 3.0);
        }
    }

    Scalar min_mass = M.minCoeff();
    return min_mass > 0;
}

} // ns mcl

#endif
