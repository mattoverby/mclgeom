// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_READELENODE_HPP
#define MCL_GEOM_READELENODE_HPP 1

#include <Eigen/Dense>
#include <fstream>
#include <vector>

namespace mcl {

static inline bool
read_ele_node(const std::string& filename_without_extensions, // will look for .node and .ele
              Eigen::MatrixXd& V,                             // vertices
              Eigen::MatrixXi& T)                             // tets
{
    using namespace Eigen;

    std::vector<Vector3d> vertices;
    std::vector<Vector4i> tets;
    const std::string& file = filename_without_extensions;

    { // Load ele

        // Load the vertices of the tetmesh
        std::stringstream ele_file;
        ele_file << file << ".ele";
        std::ifstream filestream;
        filestream.open(ele_file.str().c_str());
        if (!filestream) {
            printf("**mcl::read_elenode Error: Could not load %s\n", ele_file.str().c_str());
            return false;
        }

        std::string header;
        std::getline(filestream, header);
        std::stringstream headerSS(header);
        int n_tets = 0;
        headerSS >> n_tets;

        tets.resize(n_tets);
        std::vector<int> tet_set(n_tets, 0);
        bool starts_with_one = false;

        for (int i = 0; i < n_tets; ++i) {
            std::string line;
            std::getline(filestream, line);

            std::stringstream lineSS(line);
            size_t idx;
            int node_ids[4];
            lineSS >> idx >> node_ids[0] >> node_ids[1] >> node_ids[2] >> node_ids[3];

            // Check for 1-indexed
            if (i == 0 && idx == 1) {
                starts_with_one = true;
            }
            if (starts_with_one) {
                idx -= 1;
                for (int j = 0; j < 4; ++j) {
                    node_ids[j] -= 1;
                }
            }

            if (idx > tets.size()) {
                printf("**mcl::read_elenode Error: Indices are bad for file %s\n", ele_file.str().c_str());
                return false;
            }

            tets[idx] = Vector4i(node_ids[0], node_ids[1], node_ids[2], node_ids[3]);
            tet_set[idx] = 1;
        }
        filestream.close();

        size_t tet_set_size = tet_set.size();
        for (size_t i = 0; i < tet_set_size; ++i) {
            if (tet_set[i] == 0) {
                printf("**mcl::read_elenode Error: Indices are bad for file %s\n", ele_file.str().c_str());
                return false;
            }
        }
    }

    { // Load node

        // Load the vertices of the tetmesh
        std::stringstream node_file;
        node_file << file << ".node";
        std::ifstream filestream;
        filestream.open(node_file.str().c_str());
        if (!filestream) {
            printf("**mcl::read_elenode Error: Could not load %s\n", node_file.str().c_str());
            return false;
        }

        std::string header;
        getline(filestream, header);
        std::stringstream headerSS(header);
        int n_nodes = 0;
        headerSS >> n_nodes;

        vertices.resize(n_nodes);
        std::vector<int> vertex_set(n_nodes, 0);
        bool starts_with_one = false;

        for (int i = 0; i < n_nodes; ++i) {
            std::string line;
            std::getline(filestream, line);

            std::stringstream lineSS(line);
            double x, y, z;
            size_t idx;
            lineSS >> idx >> x >> y >> z;

            // Check for 1-indexed
            if (i == 0 && idx == 1) {
                starts_with_one = true;
            }
            if (starts_with_one) {
                idx -= 1;
            }

            if (idx > vertices.size()) {
                printf("**mcl::read_elenode Error: Indices are bad for file %s\n", node_file.str().c_str());
                return false;
            }

            vertices[idx] = Vector3d(x, y, z);
            vertex_set[idx] = 1;
        }
        filestream.close();

        size_t vert_set = vertex_set.size();
        for (size_t i = 0; i < vert_set; ++i) {
            if (vertex_set[i] == 0) {
                printf("**mcl::read_elenode Error: Indices are bad for file %s\n", node_file.str().c_str());
                return false;
            }
        }
    }

    // Check for inverted tets, and reorder if needed
    int n_tets = tets.size();
    for (int i = 0; i < n_tets; ++i) {
        Vector4i tet = tets[i];
        Vector3d a = vertices[tet[0]];
        double vv = (vertices[tet[1]] - a).dot((vertices[tet[2]] - a).cross(vertices[tet[3]] - a)) / 6.0;
        if (vv < 0) { // reorder base
            tets[i][1] = tet[2];
            tets[i][2] = tet[1];
        }
    }

    if ((int)vertices.size() == 0 || n_tets == 0) {
        printf("**mcl::read_elenode Error:  Problem loading files\n");
        return false;
    }

    // Copy to V, T
    int nv = vertices.size();
    V.resize(nv, 3);
    for (int i = 0; i < nv; ++i) {
        V.row(i) = vertices[i];
    }
    T.resize(n_tets, 4);
    for (int i = 0; i < n_tets; ++i) {
        T.row(i) = tets[i];
    }

    return true;
}

} // namespace mcl

#endif
