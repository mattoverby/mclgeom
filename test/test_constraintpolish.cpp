// Copyright Matt Overby 2021.
// Distributed under the MIT License.
#include <MCL/AssertHandler.hpp>
#include <MCL/FacesFromTets.hpp>
#include <MCL/MeshInjectivitySolver.hpp>
#include <MCL/MicroTimer.hpp>
#include <MCL/Normal.hpp>
#include <MCL/ReadEleNode.hpp>
#include <MCL/ReadVTK.hpp>
#include <MCL/SignedMeasure.hpp>
#include <MCL/Centerize.hpp>

#include <iostream>
#include <filesystem>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXd;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> RowMatrixXi;

/// @brief ConstraintZone creation and merge behavior.
void
test_constraint_zone();

/// @brief Runs FF data set
void
regression_test();

/// @brief Helper for counting inverted tetrahedra
int
count_flipped_tets(const RowMatrixXd& V, const RowMatrixXi& T,
    const std::unordered_set<int> &pinned_vertices = {})
{
    auto has_free_vertex = [&](int i) {
        for (int j = 0; j < T.cols(); ++j) {
            if (pinned_vertices.count(T(i,j) == 0)) {
                return true;
            }
        }
        return false;
    };

    int inverted_count = 0;
    for (int i = 0; i < T.rows(); ++i) {
        if (has_free_vertex(i)) {
            auto stencil = mcl::get_primitive<4>(i, T.data());
            auto v = mcl::get_verts<double, 3, 4>(V.data(), stencil.data());
            if (mcl::signed_tet_volume(v[0], v[1], v[2], v[3]) < 0.0) {
                ++inverted_count;
            }
        }
    }
    return inverted_count;
}

int
main(int argc, char* argv[])
{
    using namespace Eigen;
    (void)(argc);
    (void)(argv);

    test_constraint_zone();

    RowMatrixXd V, V0;
    RowMatrixXi T;

    // Load mesh
    {
        MatrixXd inV;
        MatrixXi inT;
        if (!mcl::read_ele_node(MCLGEOM_ROOT_DIR "/test/armadillo_3k", inV, inT)) {
            std::cout << "Failed to load " << MCLGEOM_ROOT_DIR "/test/armadillo_3k" << std::endl;
            return EXIT_FAILURE;
        }
        V0 = inV;
        V = inV;
        T = inT;
    }

    // Flip a few tets so they are inverted
    for (int i = 0; i < 10; ++i) {
        int tet_index = (i * 10) % T.rows();
        auto stencil = mcl::get_primitive<4>(tet_index, T.data());
        auto v = mcl::get_verts<double, 3, 4>(V.data(), stencil.data());
        Vector3d n = mcl::triangle_normal(v[1], v[2], v[3]);
        double h = n.dot(v[0] - v[1]);
        v[0] -= 2.0 * h * n;
        mclAssert(mcl::signed_tet_volume(v[0], v[1], v[2], v[3]) < 0.0);
        V.row(T(tet_index, 0)) = v[0];
    }

    // Solve inversions
    mcl::MeshInjectivitySolver<double, 3> solver;
    int flipped_before = count_flipped_tets(V, T);
    std::cout << "flipped before solve: " << flipped_before << std::endl;
    int iters = solver.solve(V.data(), V0.data(), V.rows(), T.data(), T.rows());
    printf("Solved in %d iterations\n", iters);

    int flipped_tets = count_flipped_tets(V, T);
    std::cout << "flipped after solve: " << flipped_tets << std::endl;
    mclAssert(flipped_tets == 0);

    regression_test();

    return EXIT_SUCCESS;
}

void
test_constraint_zone()
{
    std::vector<int> stencil_a = {0, 1, 2, 3};
    std::vector<int> stencil_b = {3, 4, 5, 6};
    std::vector<int> stencil_c = {7, 8, 9, 10};

    std::vector<mcl::ConstraintZone> zones;
    zones.emplace_back(0, stencil_a.data(), static_cast<int>(stencil_a.size()));
    zones.emplace_back(1, stencil_b.data(), static_cast<int>(stencil_b.size()));
    zones.emplace_back(2, stencil_c.data(), static_cast<int>(stencil_c.size()));

    mclAssert(zones[0].constraints.size() == 1);
    mclAssert(zones[0].global_to_local.at(0) == 0);
    mclAssert(zones[0].global_to_local.at(3) == 3);

    mcl::ConstraintZone::merge_zones(11, zones);

    mclAssert(zones.size() == 2);

    bool merged_found = false;
    bool detached_found = false;
    for (const auto& zone : zones) {
        std::unordered_set<int> zone_constraints(zone.constraints.begin(), zone.constraints.end());
        std::unordered_set<int> zone_vertices(zone.stencil.begin(), zone.stencil.end());

        if (zone_constraints.count(0) > 0 && zone_constraints.count(1) > 0) {
            merged_found = true;
            mclAssert(zone_constraints.size() == 2);
            mclAssert(zone_vertices.size() == 7);
            mclAssert(zone_vertices.count(0) > 0);
            mclAssert(zone_vertices.count(3) > 0);
            mclAssert(zone_vertices.count(6) > 0);
            mclAssert(zone.global_to_local.at(3) == 3);
        }

        if (zone_constraints.count(2) > 0) {
            detached_found = true;
            mclAssert(zone_constraints.size() == 1);
            mclAssert(zone_vertices.size() == 4);
            mclAssert(zone_vertices.count(7) > 0);
            mclAssert(zone_vertices.count(10) > 0);
        }
    }

    mclAssert(merged_found);
    mclAssert(detached_found);
}

void
regression_test()
{
    // From: https://github.com/mattoverby/mesh-data
    std::vector<std::string> surface = {
        "airplane02",  "bust02",    "centuar02", "hand02",   "laptop02", "rabbit02",    "scissor02",   "airplane03",
        "bust03",      "centuar03", "hand03",    "laptop03", "rabbit03", "scissor03",   "armadillo02", "camel02",
        "dino_skel02", "horse02",   "man02",     "santa02",  "woman02",  "armadillo03", "camel03",     "dino_skel03",
        "horse03",     "man03",     "santa03",   "woman03",  "bird02",   "cat02",       "glass02",     "lamp02",
        "octopus02",   "scape02",   "bird03",    "cat03",    "glass03",  "lamp03",      "octopus03",   "scape03"
    };

    const std::string ff_dir = MCLGEOM_ROOT_DIR "/test/mesh-data/fixedboundary/FF/3D/";

    // Check for directory:
    // https://github.com/mattoverby/mesh-data
    // Cloned to test/
    if (!std::filesystem::is_directory(ff_dir)) {
        return;
    }

    for (size_t i = 0; i < surface.size(); ++i) {
        std::string full_filename = ff_dir + "surface/" + surface[i] + "/" + surface[i] + ".vtk";
        std::string laplace_init_filename =
            ff_dir + "surface/" + surface[i] + "/initial/" + surface[i] + "_laplace.txt";
        std::string onepoint_init_filename =
            ff_dir + "surface/" + surface[i] + "/initial/" + surface[i] + "_onepoint.txt";
        std::string random_init_filename = ff_dir + "surface/" + surface[i] + "/initial/" + surface[i] + "_random.txt";

        std::cout << surface[i] << std::endl;
        std::cout << "\tloading mesh" << std::endl;

        // Load mesh
        Eigen::MatrixXd Vcm;
        Eigen::MatrixXi Tcm;
        if (!mcl::readVTK(full_filename, Vcm, Tcm)) {
            std::cerr << "Failed to load " << full_filename << std::endl;
            continue;
        }
        RowMatrixXd V0 = Vcm;
        RowMatrixXi T = Tcm;

        // Load initializers
        auto load_initializer = [&](const std::string& filename, RowMatrixXd& initializer) {
            std::ifstream init_if(filename.c_str());
            if (init_if.is_open()) {
                std::string line;
                while (std::getline(init_if, line)) {
                    std::stringstream l(line);
                    Eigen::Vector3d v(0, 0, 0);
                    int idx = -1;
                    l >> idx >> v[0] >> v[1] >> v[2];
                    if (idx >= 0 && idx < initializer.rows()) {
                        initializer.row(idx) = v;
                    }
                }
            }
        };

        std::cout << "\tloading initializers" << std::endl;

        RowMatrixXd laplace = V0;
        RowMatrixXd onepoint = V0;
        RowMatrixXd random = V0;
        load_initializer(laplace_init_filename, laplace);
        load_initializer(onepoint_init_filename, onepoint);
        load_initializer(random_init_filename, random);

        mcl::centerize(V0);
        mcl::centerize(laplace);
        mcl::centerize(onepoint);
        mcl::centerize(random);

        std::cout << "\tpinning surface" << std::endl;

        // Pin boundary verts
        Eigen::MatrixXi F;
        mcl::faces_from_tets(T, F);
        std::unordered_set<int> surface_vertices;
        for (int i = 0; i < F.rows(); ++i) {
            surface_vertices.emplace(F(i, 0));
            surface_vertices.emplace(F(i, 1));
            surface_vertices.emplace(F(i, 2));
        }
        std::vector<int> pin_inds(surface_vertices.begin(), surface_vertices.end());

        // Laplace
        if (true) {
            std::cout << "\trunning laplace initializer " << std::endl;

            int flipped_tets_init = count_flipped_tets(laplace, T, surface_vertices);
            std::cout << "\tinit flipped tets: " << flipped_tets_init << std::endl;
            mcl::MeshInjectivitySolver<double, 3> solver;
            solver.add_pins(pin_inds.data(), pin_inds.size());
            mcl::MicroTimer t;
            int iters =
                solver.solve(laplace.data(), V0.data(), V0.rows(), T.data(), T.rows());
            double ms = t.elapsed_ms();
            int flipped_tets_solved = count_flipped_tets(laplace, T, surface_vertices);
            std::cout << "\tfinal flipped tets: " << flipped_tets_solved << " in " << ms << "ms" << std::endl;
        }

        // onepoint
        if (true) {
            std::cout << "\trunning onepoint initializer " << std::endl;

            int flipped_tets_init = count_flipped_tets(onepoint, T, surface_vertices);
            std::cout << "\tinit flipped tets: " << flipped_tets_init << std::endl;
            mcl::MeshInjectivitySolver<double, 3> solver;
            solver.add_pins(pin_inds.data(), pin_inds.size());
            mcl::MicroTimer t;
            int iters =
                solver.solve(onepoint.data(), V0.data(), V0.rows(), T.data(), T.rows());
            double ms = t.elapsed_ms();
            int flipped_tets_solved = count_flipped_tets(onepoint, T, surface_vertices);
            std::cout << "\tfinal flipped tets: " << flipped_tets_solved << " in " << ms << "ms" << std::endl;
        }


        // random
        if (true) {
            std::cout << "\trunning random initializer " << std::endl;

            int flipped_tets_init = count_flipped_tets(random, T, surface_vertices);
            std::cout << "\tinit flipped tets: " << flipped_tets_init << std::endl;
            mcl::MeshInjectivitySolver<double, 3> solver;
            solver.add_pins(pin_inds.data(), pin_inds.size());
            mcl::MicroTimer t;
            int iters =
                solver.solve(random.data(), V0.data(), V0.rows(), T.data(), T.rows());
            double ms = t.elapsed_ms();
            int flipped_tets_solved = count_flipped_tets(random, T, surface_vertices);
            std::cout << "\tfinal flipped tets: " << flipped_tets_solved << " in " << ms << "ms" << std::endl;
        }
    }
}