// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_READVTK_HPP
#define MCL_READVTK_HPP 1

#include <Eigen/Core>
#include <vector>

namespace mcl {

// Reads VTK tet mesh
bool
readVTK(const std::string& vtk, Eigen::MatrixXd& V, Eigen::MatrixXi& T)
{
    using namespace Eigen;

    auto ret_err = [](const std::string& s) {
        printf("**mcl::readVTK Error: %s\n", s.c_str());
        return false;
    };

    std::vector<Vector3d> verts;
    std::vector<Vector4i> tets;

    auto file_exists = [](const std::string& n) {
        std::ifstream f(n.c_str());
        return f.good();
    };

    if (!file_exists(vtk)) {
        return ret_err("No vtk file");
    }

    // Read all lines into buffer
    std::vector<std::string> lines;
    std::ifstream vtk_file(vtk.c_str());
    if (vtk_file.is_open()) {
        std::string line;
        while (std::getline(vtk_file, line)) {
            lines.emplace_back(line);
        }
        vtk_file.close();
    }

    if (lines.size() == 0) {
        return ret_err("No data");
    }
    int num_lines = lines.size();
    int min_idx = 999;

    // Get pts and tets
    {
        std::string buff;
        int n_pts = 0;
        std::stringstream lV(lines[4]);
        lV >> buff >> n_pts;
        if (n_pts == 0) {
            return ret_err("No points");
        }
        for (int i = 5; i < n_pts + 5; ++i) {
            std::stringstream l(lines[i]);
            Vector3d v(0, 0, 0);
            l >> v[0] >> v[1] >> v[2];
            verts.emplace_back(v);
        }

        // loop until we find tets
        // then add them as needed.
        for (int i = n_pts + 5; i < num_lines; ++i) {
            std::stringstream lT(lines[i]);

            // Empty string?
            if (!lT.good())
                continue;

            // nvar should be 4
            int nvar = 0;
            lT >> nvar;
            if (nvar != 4)
                continue;

            Vector4i t(0, 0, 0, 0);
            lT >> t[0] >> t[1] >> t[2] >> t[3];
            min_idx = std::min(min_idx, t.minCoeff());
            tets.emplace_back(t);
        }
        /*
                std::stringstream lT(lines[n_pts+5]);
                lT >> buff >> n_tets;
                if (n_tets==0) { return ret_err("No tets"); }
                for (int i=n_pts+6; i<(n_pts+6)+n_tets; ++i)
                {
                    std::stringstream l(lines[i]);
                    Vector4i t(0,0,0,0);
                    int dim = 0;
                    l >> dim;
                    if( dim != 4 ){
                        return ret_err("Not a tet mesh");
                    }
                    l >> t[0] >> t[1] >> t[2] >> t[3];
                    tets.emplace_back(t);
                }
        */
    }

    if (verts.size() == 0 || tets.size() == 0) {
        return ret_err("No verts/tets");
    }

    { // Copy to V, T, start T at zero
        V.resize(verts.size(), 3);
        T.resize(tets.size(), 4);
        for (int i = 0; i < (int)verts.size(); ++i) {
            V.row(i) = verts[i];
        }
        for (int i = 0; i < (int)tets.size(); ++i) {
            T.row(i) = tets[i];
            T.row(i) -= RowVector4i::Ones() * min_idx;
        }
    }

    return verts.size() > 0 && tets.size() > 0;

} // end read vtk

} // ns mcl

#endif
