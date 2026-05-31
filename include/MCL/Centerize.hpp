// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_CENTERIZE_HPP
#define MCL_GEOM_CENTERIZE_HPP 1

#include <Eigen/Core>
#include <limits>

namespace mcl {

/// @brief Moves all of the vertices so that the center of the mesh is at the origin.
template<typename DerivedV>
inline void
centerize(Eigen::MatrixBase<DerivedV>& V);

/// @brief Returns the index of the center-most vertex
template<typename DerivedV>
inline int
get_center_index(const Eigen::MatrixBase<DerivedV>& V);

/// @brief Centerizes and scales all of the vertices in V to a target radius. Returns radius.
template<typename DerivedV>
typename DerivedV::Scalar
centerize_and_scale(Eigen::PlainObjectBase<DerivedV>& V, typename DerivedV::Scalar target_radius);

//
// Implemenation
//

template<typename DerivedV>
void
centerize(Eigen::MatrixBase<DerivedV>& V)
{
    int cols = V.cols();
    for (int i = 0; i < cols; ++i) {
        typename DerivedV::Scalar ci = V.col(i).mean();
        V.col(i).array() -= ci;
    }
} // end centerize

template<typename DerivedV>
int
get_center_index(const Eigen::MatrixBase<DerivedV>& V)
{
    typedef typename DerivedV::Scalar T;
    typedef Eigen::Matrix<T, 3, 1> Vec3t;
    int cols = std::min(3, int(V.cols()));

    Vec3t center = Vec3t::Zero();
    for (int i = 0; i < cols; ++i)
        center[i] = V.col(i).mean();

    int min_idx = -1;
    T min_dist = std::numeric_limits<T>::max();

    int nv = V.rows();
    for (int i = 0; i < nv; ++i) {
        Vec3t vi = Vec3t::Zero();
        for (int j = 0; j < cols; ++j)
            vi[j] = V(i, j);

        T dist = (center - vi).norm();
        if (dist < min_dist) {
            min_dist = dist;
            min_idx = i;
        }
    }

    return min_idx;
}

template<typename DerivedV>
typename DerivedV::Scalar
centerize_and_scale(Eigen::PlainObjectBase<DerivedV>& V, typename DerivedV::Scalar target_radius)
{
    using Scalar = typename DerivedV::Scalar;

    if (V.rows() == 0) {
        return Scalar(1);
    }

    // center at centroid
    V.rowwise() -= V.colwise().mean();

    Scalar radius = Scalar(0);
    for (Eigen::Index i = 0; i < V.rows(); ++i) {
        radius = std::max(radius, V.row(i).norm());
    }

    if (radius <= Scalar(0)) {
        return Scalar(1);
    }

    Scalar scale = target_radius / radius;
    V *= scale;

    return scale;
}

} // end ns mcl

#endif
