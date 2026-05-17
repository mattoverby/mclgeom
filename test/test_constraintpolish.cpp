// Copyright Matt Overby 2021.
// Distributed under the MIT License.
#include <MCL/ConstraintPolish.hpp>

#include <MCL/AssertHandler.hpp>
#include <MCL/MicroTimer.hpp>
#include <MCL/Normal.hpp>
#include <MCL/SignedMeasure.hpp>
#include <MCL/ReadEleNode.hpp>

#include <iostream>

int
main(int argc, char* argv[])
{
    using namespace Eigen;
    (void)(argc);
    (void)(argv);

    typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;
    typedef Matrix<int, Dynamic, Dynamic, RowMajor> RowMatrixXi;
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
    for (int i = 0; i < 10; ++i)
    {
        int tet_index = (i * 10) % T.rows();
        auto stencil = mcl::get_primitive<4>(tet_index, T.data());
        auto v = mcl::get_verts<double,3,4>(V.data(), V.data(), stencil.data(), 1.0);
        Vector3d n = mcl::triangle_normal(v[1], v[2], v[3]);
        double h = n.dot(v[0]-v[1]);
        v[0] -= 2.0 * h * n;
        mclAssert(mcl::signed_tet_volume(v[0], v[1], v[2], v[3]) < 0);
        V.row(T(tet_index,0)) = v[0];
    }

    // Gather constraints
    mcl::InequalityConstraintSet<double,3> constraint_set;
    constraint_set.add_inversions(V.data(), V0.data(), T.data(), T.rows());


    return EXIT_SUCCESS;
}
