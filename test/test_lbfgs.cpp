// Copyright Matt Overby 2021.
// Distributed under the MIT License.
#include <MCL/LBFGS.hpp>

#include <MCL/AssertHandler.hpp>
#include <MCL/Centerize.hpp>
#include <MCL/EnergyModel.hpp>
#include <MCL/FacesFromTets.hpp>
#include <MCL/ReadEleNode.hpp>
#include <MCL/SignedMeasure.hpp>

#include <Eigen/SVD>
#include <igl/edges.h>
#include <igl/opengl/glfw/Viewer.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

#include <iostream>

using namespace Eigen;

int
test_mesh_deform();

int
main(int argc, char* argv[])
{
    (void)(argc);
    (void)(argv);
    return test_mesh_deform();
}

int
test_mesh_deform()
{
    typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;
    typedef Matrix<int, Dynamic, Dynamic, RowMajor> RowMatrixXi;
    RowMatrixXd V;
    RowMatrixXi T;
    Eigen::MatrixXi F;

    // Load mesh
    {
        MatrixXd inV;
        MatrixXi inT;
        if (!mcl::read_ele_node(MCLGEOM_ROOT_DIR "/test/armadillo_3k", inV, inT)) {
            std::cout << "Failed to load " << MCLGEOM_ROOT_DIR "/test/armadillo_3k" << std::endl;
            return EXIT_FAILURE;
        }
        V = inV;
        T = inT;
        mcl::faces_from_tets(T, F);
        mcl::centerize_and_scale(V, 1);
    }

    // Test every available energy model available (except ARAP)
    std::vector<mcl::EnergyModel<double, 3>> models;
    for (int i = 0; i < mcl::ENERGY_MODEL_NUM; ++i) {
        if (i == mcl::ENERGY_MODEL_ARAP) {
            continue;
        }
        models.emplace_back(mcl::Lame<double>::soft_rubber(i));
    }
    mclAssert(!models.empty());

    // Compute rest shape
    std::vector<Eigen::Matrix<double, 3, 4>> reduction_matrices(T.rows());
    std::vector<double> volumes(T.rows(), 0);
    for (int i = 0; i < T.rows(); ++i) {
        std::array<Eigen::Vector3d, 4> xi = { V.row(T(i, 0)), V.row(T(i, 1)), V.row(T(i, 2)), V.row(T(i, 3)) };
        reduction_matrices[i] = mcl::tet_reduction_matrix(xi[0], xi[1], xi[2], xi[3]);
        volumes[i] = mcl::signed_tet_volume(xi[0], xi[1], xi[2], xi[3]);
        mclAssert(reduction_matrices[i].allFinite());
        mclAssert(volumes[i] > 0);
    }

    // Set up solver
    mcl::LBFGS<RowMatrixXd> solver;
    solver.gradient = [&](const RowMatrixXd& X, RowMatrixXd& G) -> double {
        bool compute_gradient = G.rows() == X.rows() && G.cols() == X.cols();
        if (compute_gradient) {
            G.setZero();
        }

        Eigen::VectorXd objective = Eigen::VectorXd::Zero(T.rows());
        tbb::enumerable_thread_specific<RowMatrixXd> thread_gradients(RowMatrixXd::Zero(X.rows(), 3));

        tbb::parallel_for(tbb::blocked_range<int>(0, T.rows()), [&](const tbb::blocked_range<int>& range) {
            RowMatrixXd* local_gradient = compute_gradient ? &thread_gradients.local() : nullptr;
            for (int i = range.begin(); i != range.end(); ++i) {
                const auto& energy_model = models[i % int(models.size())]; // cycle models for testing
                Eigen::Matrix<double, 4, 3> Xi;
                Xi.row(0) = X.row(T(i, 0));
                Xi.row(1) = X.row(T(i, 1));
                Xi.row(2) = X.row(T(i, 2));
                Xi.row(3) = X.row(T(i, 3));
                const Eigen::Matrix<double, 3, 4>& D = reduction_matrices[i];
                Eigen::Matrix3d F = D * Xi;

                // Hyperelastic models are only defined for positive-volume tetrahedra.
                // Returning inf will block linesearch until it finds a valid step.
                if (F.determinant() <= 0) {
                    objective[i] = std::numeric_limits<double>::infinity();
                    continue;
                }

                Eigen::JacobiSVD<Eigen::Matrix3d> svd(F, ComputeFullU | ComputeFullV);
                Eigen::Vector3d S = svd.singularValues();
                Eigen::Vector3d g = Eigen::Vector3d::Zero();
                objective[i] = energy_model.gradient(S, g) * volumes[i];

                if (compute_gradient) {
                    Eigen::Matrix3d Piola = svd.matrixU() * g.asDiagonal() * svd.matrixV().transpose();
                    Eigen::Matrix<double, 4, 3> Gi = D.transpose() * Piola * volumes[i];
                    local_gradient->row(T(i, 0)) += Gi.row(0);
                    local_gradient->row(T(i, 1)) += Gi.row(1);
                    local_gradient->row(T(i, 2)) += Gi.row(2);
                    local_gradient->row(T(i, 3)) += Gi.row(3);
                }
            }
        });

        if (compute_gradient) {
            for (const RowMatrixXd& local_gradient : thread_gradients) {
                G += local_gradient;
            }
        }

        return objective.sum();
    };

    // Compress the dillo so it springs outward.
    // Note, we are effectively doing an interior point method with returning infinity upon
    // an inverted tet. So the initial starting state must be valid - no inversions.
    for (int i = 0; i < V.rows(); ++i) {
        V(i, 1) *= 0.5;
    }

    // Init viewer
    RowMatrixXd V_prev = V;
    igl::opengl::glfw::Viewer viewer;
    viewer.core().is_animating = true;
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&) -> bool {
        V_prev = V;
        double obj = solver.iterate(V);
        std::cout << "dx: " << (V - V_prev).norm() << "\t\t energy: " << obj << std::endl;
        viewer.data().clear();
        viewer.data().set_mesh(V, F);
        return false;
    };

    viewer.launch();
    return EXIT_SUCCESS;
}