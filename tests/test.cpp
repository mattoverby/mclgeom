// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include <iostream>
#include "MCL/AssertHandler.hpp"
#include "MCL/Projection.hpp"
#include "MCL/Barycoords.hpp"
#include "MCL/LBFGS.hpp"

using namespace Eigen;

void test_projection_barys();
void test_lbfgs();

// Should probably replace with actual unit testing through CTest
int main(int argc, char *argv[])
{
    (void)(argc);
    (void)(argv);
    //test_projection_barys();
    test_lbfgs();
    return EXIT_SUCCESS;
}

void test_projection_barys()
{
    Vector3d f[3] = {
        Vector3d::Zero(),
        Vector3d(1,0,0),
        Vector3d(0,1,0) };
    Vector3d pt(2,2,2);
    Vector3d proj = mcl::point_on_triangle(pt, f[0], f[1], f[2]);
    Vector3d barys = mcl::point_triangle_barys(pt, f[0], f[1], f[2]);
    Vector3d proj_barys = mcl::point_triangle_barys(proj, f[0], f[1], f[2]);
    mclAssert(barys.maxCoeff() > 1);
    mclAssert(proj_barys.maxCoeff() < 1);
    std::cout << barys.transpose() << " " << proj_barys.transpose() << std::endl;
}

void test_lbfgs()
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

}
