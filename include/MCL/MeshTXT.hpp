// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef GINI_MESHTXT_HPP
#define GINI_MESHTXT_HPP 1

#include <Eigen/Geometry>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace mcl {

// Simple, slow, plain text
//
// X is n x DIM vertices (DIM = 2 or 3)
// P is m x PDIM primitives (PDIM = 3 for tris, 4 for tets)
// Returns true on success
// TODO: some error checking on read

static inline bool
write_mesh_txt(const std::string& filename, const Eigen::MatrixXd& X, const Eigen::MatrixXi& P);

static inline bool
read_mesh_txt(const std::string& filename, const Eigen::MatrixXd& X, const Eigen::MatrixXi& P);

//
// Implementation
//

bool
write_mesh_txt(const std::string& filename, const Eigen::MatrixXd& X, const Eigen::MatrixXi& P)
{
    char delim = ' ';

    if (X.cols() < 2 || X.cols() > 3) {
        printf("Cannot write mesh txt, bad dim X\n");
        return false;
    }

    if (P.cols() < 3 || P.cols() > 4) {
        printf("Cannot write mesh txt, bad dim P\n");
        return false;
    }

    std::ofstream fi(filename);
    fi << std::fixed << std::setprecision(16);

    int nx = X.rows();
    int xc = X.cols();
    for (int i = 0; i < nx; ++i) {
        fi << 'v';
        for (int j = 0; j < xc; ++j) {
            fi << delim << X(i, j);
        }
        fi << std::endl;
    }

    int np = P.rows();
    int pc = P.cols();
    char prim_tag = pc == 3 ? 'f' : 't'; // face, tet
    for (int i = 0; i < np; ++i) {
        fi << prim_tag;
        for (int j = 0; j < pc; ++j) {
            fi << delim << P(i, j);
        }
        fi << std::endl;
    }

    fi.close();
    return true;
}

bool
read_mesh_txt(const std::string& filename, Eigen::MatrixXd& X, Eigen::MatrixXi& P)
{
    int x_dim = 0;
    int p_dim = 0;
    using namespace Eigen;
    std::vector<Vector3d> x;
    std::vector<Vector4i> p;

    std::ifstream fi(filename);
    if (fi.is_open()) {
        std::string line;
        while (std::getline(fi, line)) {
            std::stringstream ss(line);
            char tag = 'X';
            ss >> tag;
            switch (tag) {
                default: {
                    printf("Bad line: %s\n", line.c_str());
                } break;
                case 'v': {
                    x.emplace_back(Vector3d::Zero());
                    for (int i = 0; i < 3 && ss.good(); ++i) {
                        x_dim = std::max(x_dim, i + 1);
                        ss >> x.back()[i];
                    }
                } break;
                case 'f': {
                    p.emplace_back(Vector4i::Zero());
                    for (int i = 0; i < 3 && ss.good(); ++i) {
                        p_dim = std::max(p_dim, i + 1);
                        ss >> p.back()[i];
                    }
                } break;
                case 't': {
                    p.emplace_back(Vector4i::Zero());
                    for (int i = 0; i < 4 && ss.good(); ++i) {
                        p_dim = std::max(p_dim, i + 1);
                        ss >> p.back()[i];
                    }
                } break;
            }
        }
        fi.close();
    } else {
        printf("Could not open %s\n", filename.c_str());
        return false;
    }

    // Make X
    int nx = x.size();
    X.resize(nx, x_dim);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < x_dim; ++j) {
            X(i, j) = x[i][j];
        }
    }

    // Make P
    int np = p.size();
    P.resize(np, p_dim);
    for (int i = 0; i < np; ++i) {
        for (int j = 0; j < p_dim; ++j) {
            P(i, j) = p[i][j];
        }
    }

    return true;
}

} // namespace mcl

#endif