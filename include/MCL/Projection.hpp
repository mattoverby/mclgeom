// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_PROJECTION_HPP
#define MCL_PROJECTION_HPP 1

#include <Eigen/Core>

namespace mcl
{


// Projection on Triangle
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_triangle(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &p1, const Eigen::Matrix<T,3,1> &p2, const Eigen::Matrix<T,3,1> &p3);

// Projection on Sphere
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_sphere(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &center, const T &rad);

// Projection on a Box
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_box(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &bmin, const Eigen::Matrix<T,3,1> &bmax);

// Project a point on to a plane
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_plane(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &plane_norm, const Eigen::Matrix<T,3,1> &plane_pt);

// Projection on an edge
template <typename T> 
static inline Eigen::Matrix<T,2,1> point_on_edge(const Eigen::Matrix<T,2,1> &point, const Eigen::Matrix<T,2,1> &p1, const Eigen::Matrix<T,2,1> &p2);

// Projection on an edge (3D)
template <typename T>
static inline Eigen::Matrix<T,3,1> point_on_edge(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &p1, const Eigen::Matrix<T,3,1> &p2);

//
//	Implementation
//

// I do not know the source of the function but have seen it around in various libs
template <typename T>
Eigen::Matrix<T,3,1> point_on_triangle(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &p1, const Eigen::Matrix<T,3,1> &p2, const Eigen::Matrix<T,3,1> &p3)
{
	auto clamp_zero_one = [](const T &val){ return val < 0 ? 0 : (val > 1 ? 1 : val); };

	Eigen::Matrix<T,3,1> edge0 = p2 - p1;
	Eigen::Matrix<T,3,1> edge1 = p3 - p1;
	Eigen::Matrix<T,3,1> v0 = p1 - point;

	T a = edge0.dot(edge0);
	T b = edge0.dot(edge1);
	T c = edge1.dot(edge1);
	T d = edge0.dot(v0);
	T e = edge1.dot(v0);
	T det = a*c - b*b;
	T s = b*e - c*d;
	T t = b*d - a*e;

	const T zero(0);
	const T one(1);

	if ( s + t < det ) {
		if ( s < zero ) {
		    if ( t < zero ) {
			if ( d < zero ) {
			    s = clamp_zero_one( -d/a );
			    t = zero;
			}
			else {
			    s = zero;
			    t = clamp_zero_one( -e/c );
			}
		    }
		    else {
			s = zero;
			t = clamp_zero_one( -e/c );
		    }
		}
		else if ( t < zero ) {
		    s = clamp_zero_one( -d/a );
		    t = zero;
		}
		else {
		    T invDet = one / det;
		    s *= invDet;
		    t *= invDet;
		}
	}
	else {
		if ( s < zero ) {
		    T tmp0 = b+d;
		    T tmp1 = c+e;
		    if ( tmp1 > tmp0 ) {
			T numer = tmp1 - tmp0;
			T denom = a-T(2)*b+c;
			s = clamp_zero_one( numer/denom );
			t = one-s;
		    }
		    else {
			t = clamp_zero_one( -e/c );
			s = zero;
		    }
		}
		else if ( t < zero ) {
		    if ( a+d > b+e ) {
			T numer = c+e-b-d;
			T denom = a-T(2)*b+c;
			s = clamp_zero_one( numer/denom );
			t = one-s;
		    }
		    else {
			s = clamp_zero_one( -e/c );
			t = zero;
		    }
		}
		else {
		    T numer = c+e-b-d;
		    T denom = a-T(2)*b+c;
		    s = clamp_zero_one( numer/denom );
		    t = one - s;
		}
	}

	return (p1 + edge0*s + edge1*t);

} // end project triangle


template <typename T>
Eigen::Matrix<T,3,1> point_on_sphere(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &center, const T &rad)
{
	Eigen::Matrix<T,3,1> dir = point-center;
	dir.normalize();
	return (center + dir*rad);
} // end project sphere


template <typename T>
Eigen::Matrix<T,3,1> point_on_box(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &bmin, const Eigen::Matrix<T,3,1> &bmax)
{
	// Loops through axes and moves point to nearest surface
	Eigen::Matrix<T,3,1> x = point;
	T dx = std::numeric_limits<T>::max();
	for (int i=0; i<3; ++i)
	{
		T dx_max = std::abs(bmax[i]-point[i]);
		T dx_min = std::abs(bmin[i]-point[i]);
		if(dx_max < dx)
		{
			x = point;
			x[i] = bmax[i];
			dx = dx_max;
		}
		if (dx_min < dx)
		{
			x = point;
			x[i] = bmin[i];
			dx = dx_min;
		}
	}
	return x;
} // end project box


template <typename T>
Eigen::Matrix<T,3,1> point_on_plane(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &plane_norm, const Eigen::Matrix<T,3,1> &plane_pt)
{
    T d = -1 * plane_norm.dot(point-plane_pt);
    Eigen::Matrix<T,3,1> t_vec = plane_norm * d;
    return point + t_vec;
}


template <typename T>
static Eigen::Matrix<T,2,1> point_on_edge(const Eigen::Matrix<T,2,1> &p, const Eigen::Matrix<T,2,1> &e0, const Eigen::Matrix<T,2,1> &e1)
{
	Eigen::Matrix<T,2,1> e = (e1-e0);
	T e_len2 = e.dot(e);
	if(e_len2 <= 0.0) { return e0; } // zero length edge
	Eigen::Matrix<T,2,1> pe0 = (p-e0);
	T t = pe0.dot(e)/e_len2;
	if (t < 0.0){ return e0; }
	else if (t > 1.0) { return e1; }
	return e0 + t * e;
}


// Ericson, Real-Time Collision Detection
template <typename T>
static Eigen::Matrix<T,3,1> point_on_edge(const Eigen::Matrix<T,3,1> &point, const Eigen::Matrix<T,3,1> &p1, const Eigen::Matrix<T,3,1> &p2)
{
	Eigen::Matrix<T,3,1> ab = p2-p1;
	double t = ab.dot(point-p1);
	Eigen::Matrix<T,3,1> d = point; // result

	// c projects outside the [a,b] interval, on the a side; clamp to a
	if (t <= 0)
	{
		t = 0;
		d = p1;
	}
	else
	{
		T denom = ab.dot(ab); // Always nonnegative since denom = ||ab|| ∧ 2
		// c projects outside the [a,b] interval, on the b side; clamp to b
		if (t >= denom)
		{
			t = 1;
			d = p2;
		}
		// c projects inside the [a,b] interval; must do deferred divide now
		else
		{
			t = t/denom;
			d = p1+t*ab;
		}
	}
	return d;
}

} // end namespace mcl

#endif

