// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#include <iostream>
#include "MCL/AssertHandler.hpp"
#include "MCL/Projection.hpp"
#include "MCL/Barycoords.hpp"
#include "MCL/ArgParser.hpp"

using namespace Eigen;

void test_projection_barys();
void test_argparser(int argc, char *argv[]);

// Should probably replace with actual unit testing through CTest
int main(int argc, char *argv[])
{
    test_projection_barys();
    test_argparser(argc, argv);
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
    
    Vector3d p0(0,0,0);
    Vector3d p1(1,0,0);
    Vector3d q0(0,1,0);
    Vector3d q1(1,1,0);
    Vector4d ee_bary(0,0,0,0);
    Vector3d line = mcl::edge_to_edge<double>(p0, p1, q0, q1, ee_bary);
    std::cout << "Line: " << line.transpose() << ", barys: " << ee_bary.transpose() << std::endl;
}

void test_argparser(int argc, char *argv[])
{
    mcl::ArgParser args(argc, argv);
    args.save_to_file("args.txt");
}
