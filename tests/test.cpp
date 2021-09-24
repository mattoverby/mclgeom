// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include <iostream>
#include "MCL/AssertHandler.hpp"
#include "MCL/Projection.hpp"
#include "MCL/Barycoords.hpp"

using namespace Eigen;

void test_projection_barys();

// Should probably replace with actual unit testing through CTest
int main(int argc, char *argv[])
{
    (void)(argc);
    (void)(argv);
    test_projection_barys();
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