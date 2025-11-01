// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include <iostream>

#include <MCL/AssertHandler.hpp>
#include <MCL/LBFGS.hpp>
#include <MCL/ReadEleNode.hpp>
#include <MCL/Centerize.hpp>

#include <igl/edges.h>
#include <igl/opengl/glfw/Viewer.h>

using namespace Eigen;

int test_mesh_deform();

int main(int argc, char *argv[])
{
    (void)(argc);
    (void)(argv);
    return test_mesh_deform();
}

int test_mesh_deform()
{
    typedef Matrix<double,Dynamic,Dynamic,RowMajor> RowMatrixXd;
    typedef Matrix<int,Dynamic,Dynamic,RowMajor> RowMatrixXi;
    RowMatrixXd V0, V1;
    RowMatrixXi T, E;

    // Load mesh
    {
        MatrixXd inV;
        MatrixXi inT;
        if (!mcl::read_ele_node(MCLGEOM_ROOT_DIR "/test/armadillo_3k", inV, inT))
        {
            std::cout << "Failed to load " << MCLGEOM_ROOT_DIR "/test/armadillo_3k" << std::endl;
            return EXIT_FAILURE;
        }

        V0 = inV;
        V1 = inV;
        T = inT;
    }

    igl::edges(T, E);
    int ne = E.rows();
    mcl::centerize(V0);
    mcl::scale_to_sphere(V0, 1);
    V1 = V0;
    for (int i=0; i<V1.rows(); ++i)
        V1.row(i) += RowVector3d::Ones() * (double(i)/double(V1.rows()));

    // Initialize solver
    mcl::LBFGS<RowMatrixXd> solver;
    solver.gradient = [&](const RowMatrixXd& x, RowMatrixXd& g)->double
    {
        if (g.rows()>0)
            g.setZero();

        double tot_energy = 0;
        double k = 100;
        for (int i=0; i<ne; ++i)
        {
            int e0 = E(i,0);
            int e1 = E(i,1);
            RowVector3d edge = (x.row(e0)-x.row(e1));
            double l = edge.norm();
            double r = (V0.row(e0)-V0.row(e1)).norm(); // rest
            if (g.rows() == x.rows() && std::abs(l) > 1e-12)
            {
                edge /= l;
                edge *= (l-r)*k;
                mclAssert(edge.allFinite());
                g.row(e0) += edge;
                g.row(e1) -= edge;
            }
            double energy = 0.5 * k * (l-r) * (l-r); // (k/2)||l-l0||^2
            tot_energy += energy;
            mclAssert(std::isfinite(energy));
        }
        return tot_energy;
    };
    solver.initialize(V1);

    // Function to update render lines
    MatrixXd E0(ne,3), E1(ne,3), EC(ne,3);
    auto update_edge_points = [&](RowMatrixXd &V)
    {
        for (int i=0; i<ne; ++i)
        {
            E0.row(i) = V.row(E(i,0));
            E1.row(i) = V.row(E(i,1));
            EC.row(i) = RowVector3d(1,0,0);
        }
    };

    // Init viewer
    igl::opengl::glfw::Viewer viewer;
    update_edge_points(V1);
    viewer.data().add_edges(E0, E1, EC);
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool
    {
        RowMatrixXd V_prev = V1;
        double obj = solver.iterate(V1);
        std::cout << "dx: " << (V1 - V_prev).norm() << "\t\t energy: " << obj << std::endl;
        update_edge_points(V1);
        viewer.data().clear();
        viewer.data().add_edges(E0, E1, EC);
        return false;
    };

    viewer.launch();
    return EXIT_SUCCESS;
}