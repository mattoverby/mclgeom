// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include <iostream>
#include "MCL/AssertHandler.hpp"
#include "MCL/LBFGS.hpp"
#include <igl/readMSH.h>
#include <igl/edges.h>
#include <igl/opengl/glfw/Viewer.h>
#include "MCL/Centerize.hpp"

using namespace Eigen;

int test_lbfgs();
int test_beam();

// Should probably replace with actual unit testing through CTest
int main(int argc, char *argv[])
{
    (void)(argc);
    (void)(argv);
    //return test_lbfgs();
    return test_beam();
}

int test_lbfgs()
{
    mcl::LBFGS<Vector2d> solver;
    solver.gradient = [&](const Eigen::Vector2d& x, Eigen::Vector2d& g)->double
    {
        if (g.rows() == x.rows())
        {
            g[0] = 400*x[0]*x[0]*x[0] - 400*x[0]*x[1] + 2*x[0] - 2;
            g[1] = 200*(x[1] - x[0]*x[0]);
        }
        double a = 1.0 - x[0];
        double b = x[1] - x[0]*x[0];
        return a*a + b*b*100.0;
    };
    
    Vector2d x = Vector2d::Random();
    std::cout << "x init: " << x.transpose() << std::endl;
    double obj = solver.minimize(x);
    std::cout << "obj: " << obj << ", x: " << x.transpose() << std::endl;
    return EXIT_SUCCESS;
}

int test_beam()
{
    typedef Matrix<double,Dynamic,Dynamic,RowMajor> RowMatrixXd;
    typedef Matrix<int,Dynamic,Dynamic,RowMajor> RowMatrixXi;
    RowMatrixXd V0, V1;
    RowMatrixXi T, E;

    // Load mesh
    {
        MatrixXd inV;
        MatrixXi inT, inF;
        VectorXi inTT, inFT;
        if (!igl::readMSH(MCL_GEOM_ROOT_DIR "/tests/bunny32k.msh", inV, inF, inT, inFT, inTT))
        {
            std::cout << "Failed to load " << MCL_GEOM_ROOT_DIR "/tests/bunny32k.msh" << std::endl;
            return EXIT_FAILURE;
        }

        V0 = inV;
        V1 = inV;
        T = inT;
    }

    igl::edges(T, E);
    mcl::centerize(V0);
    V1 = V0;

    // Initialize solver
    mcl::LBFGS<RowMatrixXd> solver;
    solver.gradient = [&](const RowMatrixXd& x, RowMatrixXd& g)->double
    {
        return 0;
    };
    solver.initialize(V1);

    // Function to update render lines
    int ne = E.rows();
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
        solver.iterate(V1);
        update_edge_points(V1);
        viewer.data().clear();
        viewer.data().add_edges(E0, E1, EC);
        return false;
    };

    std::cout << "Press A to start solver" << std::endl;
    viewer.launch();
    return EXIT_SUCCESS;
}