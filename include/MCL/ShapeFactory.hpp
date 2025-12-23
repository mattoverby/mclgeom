// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_SHAPEFACTORY_HPP
#define MCL_GEOM_SHAPEFACTORY_HPP

#include <Eigen/Dense>

namespace mcl {

/// @brief Makes a 3D box between lower and upper corners
template<typename T, typename DerivedV, typename DerivedF>
static inline void
make_tri_box(const Eigen::Vector3<T>& bmin, const Eigen::Vector3<T>& bmax, DerivedV& V, DerivedF& F);

/// @brief Makes a 3D quad between lower and upper corners
template<typename T, typename DerivedV, typename DerivedF>
static inline void
make_quad(const Eigen::Vector2<T>& bottom_left,
          const Eigen::Vector2<T>& upper_right,
          int tessellation,
          DerivedV& V,
          DerivedF& F);

//
// Implementation
//

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

template<typename T, typename DerivedV, typename DerivedF>
void
make_quad(const Eigen::Vector2<T>& bottom_left,
          const Eigen::Vector2<T>& upper_right,
          int tessellation,
          DerivedV& V,
          DerivedF& F)
{
    V.resize((tessellation + 1) * (tessellation + 1), 3);
    F.resize(tessellation * tessellation * 2, 3);
    Eigen::Vector2<T> size = upper_right - bottom_left;
    Eigen::Vector2<T> step = size / static_cast<T>(tessellation);

    // Generate vertices
    int vertex_count = 0;
    for (int j = 0; j <= tessellation; ++j) {
        for (int i = 0; i <= tessellation; ++i) {
            Eigen::Vector2<T> p2 = bottom_left + Eigen::Vector2<T>(i * step[0], j * step[1]);
            V.row(vertex_count++) = Eigen::RowVector3<T>(p2[0], p2[1], 0);
        }
    }

    auto idx = [tessellation](int i, int j) { return j * (tessellation + 1) + i; };

    // Generate faces (two triangles per quad)
    int face_count = 0;
    for (int j = 0; j < tessellation; ++j) {
        for (int i = 0; i < tessellation; ++i) {
            int i0 = idx(i, j);
            int i1 = idx(i + 1, j);
            int i2 = idx(i, j + 1);
            int i3 = idx(i + 1, j + 1);
            F.row(face_count++) = Eigen::RowVector3i(i0, i1, i3);
            F.row(face_count++) = Eigen::RowVector3i(i0, i3, i2);
        }
    }
}

} // end namespace mcl

#endif
