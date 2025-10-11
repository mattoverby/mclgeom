// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_READSEAMEDOBJ_HPP
#define MCL_READSEAMEDOBJ_HPP 1

#include <Eigen/Core>
#include <vector>

// TODO my own readOBJ to remove igl dependency
#include <igl/readOBJ.h>

namespace mcl {

static inline bool
read_seamed_obj(std::string filename,
                Eigen::MatrixXd& V3D, // 3D vertices
                Eigen::MatrixXd& VTC, // 2D UV tex init, if found
                Eigen::MatrixXi& F)
{
    using namespace Eigen;
    MatrixXd V3D_, vN_unused;
    MatrixXi F3D_, fN_unused;
    bool read_success = igl::readOBJ(filename, V3D_, VTC, vN_unused, F3D_, F, fN_unused);

    if (!read_success) {
        return false;
    }

    if (VTC.rows() > 0 && F.rows() > 0) {
        // Create 3D verts from F so face meshes match
        std::vector<bool> filled(VTC.rows(), false);
        V3D = MatrixXd::Zero(VTC.rows(), 3);
        if (F3D_.rows() != F.rows()) {
            return false;
        }
        int nf = F3D_.rows();
        for (int i = 0; i < nf; i++) {
            for (int j : { 0, 1, 2 }) {
                int f = F3D_(i, j);
                int f_uv = F(i, j);
                if (!filled[f_uv]) {
                    V3D.row(f_uv) = V3D_.row(f);
                    filled[f_uv] = true;
                }
            }
        }
    } else {
        F = F3D_;
        VTC = V3D_.block(0, 0, V3D_.rows(), 2);
        V3D = V3D_;
    }

    if (F.minCoeff() > 0) {
        F.array() -= 1;
    }
    return true;
}

} // end ns mcl

#endif