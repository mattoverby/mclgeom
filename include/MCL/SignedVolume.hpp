// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_SIGNEDVOLUME_HPP
#define MCL_SIGNEDVOLUME_HPP 1

#include <Eigen/Dense>

namespace mcl
{

static inline double signed_triangle_area(
	const Eigen::Vector2d &p1,
	const Eigen::Vector2d &p2,
	const Eigen::Vector2d &p3);

// Returns (unscaled) first derivative of signed triangle area
static inline std::vector<Eigen::Vector2d> signed_triangle_area_gradients(
	const Eigen::Vector2d &p1,
	const Eigen::Vector2d &p2,
	const Eigen::Vector2d &p3);

static inline double triangle_area(
	const Eigen::Vector3d &p1,
	const Eigen::Vector3d &p2,
	const Eigen::Vector3d &p3);

static inline double signed_tet_volume(
	const Eigen::Vector3d &p1,
	const Eigen::Vector3d &p2,
	const Eigen::Vector3d &p3,
	const Eigen::Vector3d &p4);

// Returns (unscaled) first derivative of signed tet volume
static inline std::vector<Eigen::Vector3d> signed_tet_volume_gradients(
	const Eigen::Vector3d &p1,
	const Eigen::Vector3d &p2,
	const Eigen::Vector3d &p3,
	const Eigen::Vector3d &p4);

static inline double tet_surface_area(
	const Eigen::Vector3d &p1,
	const Eigen::Vector3d &p2,
	const Eigen::Vector3d &p3,
	const Eigen::Vector3d &p4);

// Probably not where this belongs but oh well.
// Returns the faces of a tet
static inline std::vector<Eigen::Vector3i>
	faces_from_tet(const Eigen::RowVector4i &t);

//
// Implementation
//

inline double signed_triangle_area(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Vector2d &p3)
{
	return 0.5 * ( -p2[0]*p1[1] + p3[0]*p1[1] + p1[0]*p2[1] - p3[0]*p2[1] - p1[0]*p3[1] + p2[0]*p3[1] );
}

inline std::vector<Eigen::Vector2d> signed_triangle_area_gradients(
	const Eigen::Vector2d &a, const Eigen::Vector2d &b, const Eigen::Vector2d &c)
{
	std::vector<Eigen::Vector2d> g(3);
	g[0] = 0.5 * Eigen::Vector2d(b[1]-c[1], -b[0]+c[0]);
	g[1] = 0.5 * Eigen::Vector2d(-a[1]+c[1], a[0]-c[0]);
	g[2] = 0.5 * Eigen::Vector2d(a[1]-b[1], -a[0]+b[0]);
	return g;
}

// https://en.wikipedia.org/wiki/Heron%27s_formula
inline double triangle_area(const Eigen::Vector3d &p1, const Eigen::Vector3d &p2, const Eigen::Vector3d &p3)
{
	double a = (p1-p2).norm();
	double b = (p2-p3).norm();
	double c = (p3-p1).norm();
	double s = (a+b+c) * 0.5;
	return std::sqrt(s*(s-a)*(s-b)*(s-c));
}

inline double signed_tet_volume(
	const Eigen::Vector3d &p1, const Eigen::Vector3d &p2,
	const Eigen::Vector3d &p3, const Eigen::Vector3d &p4)
{
	Eigen::Matrix3d edges;
	edges.col(0) = p2 - p1;
	edges.col(1) = p3 - p1;
	edges.col(2) = p4 - p1;
	return (1.0/6.0) * edges.determinant();
}

inline std::vector<Eigen::Vector3d> signed_tet_volume_gradients(
	const Eigen::Vector3d &a, const Eigen::Vector3d &b,
	const Eigen::Vector3d &c, const Eigen::Vector3d &d)
{
	std::vector<Eigen::Vector3d> grads(4);
	const Eigen::Vector3d &p0 = a;
	const Eigen::Vector3d &p1 = b;
	const Eigen::Vector3d &p2 = c;
	const Eigen::Vector3d &p3 = d;
	static const double sixth = (1.0/6.0);
	grads[0] = sixth * (p1 - p2).cross(p3 - p2);
	grads[1] = sixth * (p2 - p0).cross(p3 - p0);
	grads[2] = sixth * (p0 - p1).cross(p3 - p1);
	grads[3] = sixth * (p1 - p0).cross(p2 - p0);
	return grads;
}

inline double tet_surface_area(
		const Eigen::Vector3d &p1,
		const Eigen::Vector3d &p2,
		const Eigen::Vector3d &p3,
		const Eigen::Vector3d &p4)
{
	double a1 = triangle_area(p2,p3,p4);
	double a2 = triangle_area(p2,p3,p1);
	double a3 = triangle_area(p3,p4,p1);
	double a4 = triangle_area(p4,p2,p1);
	return (a1+a2+a3+a4);
}

inline std::vector<Eigen::Vector3i> faces_from_tet(const Eigen::RowVector4i &t)
{
	using namespace Eigen;
	std::vector<Vector3i> f = {
		Vector3i(t[0], t[1], t[3]),
		Vector3i(t[0], t[2], t[1]),
		Vector3i(t[0], t[3], t[2]),
		Vector3i(t[1], t[2], t[3]) };
	return f;
}

} // ns mcl

#endif
