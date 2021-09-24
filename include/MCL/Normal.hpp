// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_NORMAL_HPP
#define MCL_NORMAL_HPP 1

#include <Eigen/Dense>

namespace mcl
{

template <typename T>
static inline Eigen::Matrix<T,3,1> triangle_normal(
	const Eigen::Matrix<T,3,1> &a,
	const Eigen::Matrix<T,3,1> &b,
	const Eigen::Matrix<T,3,1> &c,
	bool normalize = true)
{
	Eigen::Matrix<T,3,1> n = (b-a).cross(c-a);
	if (normalize) { n.stableNormalize(); }
	return n;
}

template <typename T>
static inline Eigen::Matrix<T,2,1> edge_normal(
	const Eigen::Matrix<T,2,1> &p0,
	const Eigen::Matrix<T,2,1> &p1,
	bool normalize = true)
{
	Eigen::Matrix<T,2,1> n(p1[1]-p0[1], -(p1[0]-p0[0]));
	if (normalize) { n.stableNormalize(); }
	return n;
}

} // ns mcl

#endif
