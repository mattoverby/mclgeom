// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_SHAPEFACTORY_HPP
#define MCL_GEOM_SHAPEFACTORY_HPP

#include <Eigen/Dense>

namespace mcl {

template<typename T, typename DerivedV, typename DerivedF>
static inline void
make_tri_box(const Eigen::Vector3<T>& bmin, const Eigen::Vector3<T>& bmax, DerivedV& V, DerivedF& F)
{
    V.resize(8, 3);
    V << bmin[0], bmin[1], bmin[2], // 0
        bmax[0], bmin[1], bmin[2],  // 1
        bmax[0], bmax[1], bmin[2],  // 2
        bmin[0], bmax[1], bmin[2],  // 3
        bmin[0], bmin[1], bmax[2],  // 4
        bmax[0], bmin[1], bmax[2],  // 5
        bmax[0], bmax[1], bmax[2],  // 6
        bmin[0], bmax[1], bmax[2];  // 7
    F.resize(12, 3);
    F << 0, 2, 1, 0, 3, 2, 4, 5, 6, 4, 6, 7, 0, 7, 3, 0, 4, 7, 1, 2, 6, 1, 6, 5, 3, 7, 6, 3, 6, 2, 0, 1, 5, 0, 5, 4;
}

} // end namespace mcl

#endif
