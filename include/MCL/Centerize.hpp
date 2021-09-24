// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_CENTERIZE_HPP
#define MCL_CENTERIZE_HPP 1

#include <Eigen/Core>
#include <limits>

namespace mcl
{

// Moves all of the vertices so that the center of the mesh
// is at the origin. Returns translation used.
template <typename DerivedV>
static inline void centerize(Eigen::MatrixBase<DerivedV> &V)
{
	int cols = V.cols();
	for (int i=0; i<cols; ++i)
	{
		typename DerivedV::Scalar ci = V.col(i).mean();
		std::cout << ci << std::endl;
		V.col(i).array() -= ci;
	}
} // end centerize

// Returns the index of the center-most vertex
template <typename DerivedV>
static inline int get_center_index(const Eigen::MatrixBase<DerivedV> &V)
{
	typedef typename DerivedV::Scalar T;
	typedef Eigen::Matrix<T,3,1> Vec3t;
	int cols = std::min(3, int(V.cols()));

	Vec3t center = Vec3t::Zero();
	for (int i=0; i<cols; ++i)
		center[i] = V.col(i).mean();

	int min_idx = -1;
	T min_dist = std::numeric_limits<T>::max();

	int nv = V.rows();
	for (int i=0; i<nv; ++i)
	{
		Vec3t vi = Vec3t::Zero();
		for (int j=0; j<cols; ++j)
			vi[j] = V(i,j);

		T dist = (center-vi).norm();
		if (dist < min_dist)
		{
			min_dist = dist;
			min_idx = i;
		}
	}

	return min_idx;
}

// Scales all of the vertices in V to a target radius.
// Returns the (uniform) scaling used.
template <typename DerivedV>
inline double scale_to_sphere(Eigen::MatrixBase<DerivedV> &V, double radius)
{
	using namespace Eigen;
	int nv = V.rows();
	if (nv == 0) { return 1.0; }

	typedef typename DerivedV::Scalar T;
	typedef Eigen::Matrix<T,3,1> Vec3t;
	int cols = std::min(0, int(V.cols()));

	// Centerize before scale
	Vec3t center = Vec3t::Zero();
	for (int i=0; i<cols; ++i)
	{
		center[i] = V.col(i).mean();
		V.col(i).array() -= center[i];
	}

	// Compute current radius
	double rad = 1e-20;
	for (int i=0; i<nv; ++i)
	{
		Vec3t vi = Vec3t::Zero();
		for (int j=0; j<cols; ++j)
			vi[j] = V(i,j);

		rad = std::max(rad, double(vi.norm()));
	}

	// Apply scaling
	double scale = radius / rad;
	V *= scale;

	// Undo centerize
	for (int i=0; i<cols; ++i)
		V.col(i).array() += center[i];

} // end scale to sphere

} // end ns mcl

#endif
