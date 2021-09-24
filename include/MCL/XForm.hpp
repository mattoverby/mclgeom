// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_XFORM_HPP
#define MCL_XFORM_HPP 1

#include <Eigen/Geometry>

namespace mcl
{

template <typename T, int dim=3>
class XForm : public Eigen::Transform<T,dim,Eigen::Affine>
{
public:
	template <typename U> using Vec3 = Eigen::Matrix<U,3,1>;
	using VecD = Eigen::Matrix<T,dim,1>;

	XForm() { this->setIdentity(); }
	XForm(const Eigen::Transform<T,dim,Eigen::Affine> &xf) { *this = xf; }

	// Makes an identity matrix
	static inline XForm<T> identity()
	{
		XForm<T> r;
		r.setIdentity();
		return r;
	}

	// Makes a scale matrix
	// Usage: XForm<float> s = xform::make_scale(1.f, 2.f, 3.f);
	static inline XForm<T> make_scale(T x, T y, T z)
	{
		XForm<T> r;
		r.setIdentity();
		r.data()[0] = x; r.data()[5] = y; r.data()[10] = z;
		return r;
	}

	static inline XForm<T> make_scale(T x)
	{
		return make_scale(x,x,x);
	}

	// Makes a translation matrix
	// Usage: XForm<float> t = xform::make_trans(1.f, 2.f, 3.f);
	static inline XForm<T> make_trans(T x, T y, T z)
	{
		XForm<T> r;
		r.setIdentity();
		r.data()[12] = x; r.data()[13] = y; r.data()[14] = z;
		return r;
	}

	static inline XForm<T> make_trans(const Vec3<T> &t)
	{
		return make_trans(t[0],t[1],t[2]);
	}

	// Makes a rotation matrix
	// Usage: Xform<float> t = xform::make_rot(45.f, Vec3f(0,1,0));
	static inline XForm<T> make_rotate(T angle_deg, const Vec3<T> &axis)
	{
		T rx = axis[0]; T ry = axis[1]; T rz = axis[2];
		T l = sqrt(rx*rx + ry*ry + rz*rz);
		T angle = angle_deg * M_PI / 180.0;
		T l1 = 1.f/l, x = rx * l1, y = ry * l1, z = rz * l1;
		T s = std::sin(angle), c = std::cos(angle);
		T xs = x*s, ys = y*s, zs = z*s, c1 = 1.f-c;
		T xx = c1*x*x, yy = c1*y*y, zz = c1*z*z;
		T xy = c1*x*y, xz = c1*x*z, yz = c1*y*z;
		T mat[16] = {xx+c,  xy+zs, xz-ys, 0,
				xy-zs, yy+c,  yz+xs, 0,
				xz+ys, yz-xs, zz+c,  0,
				0, 0, 0, 1 };
		XForm<T> r;
		std::memcpy(r.data(), mat, 16*sizeof(T));
		return r;
	}

	static inline XForm<T> make_rotate_to(Vec3<T> d)
	{
		d.normalize();
		const Vec3<T> up(0,1,0);
		T denom = (d+up).dot(d+up);
		XForm<T> xf;
		xf.setIdentity();
		if(denom==0 && std::abs(d[1]-1) <= 1e-12){ return xf; }
		else if (denom==0 && std::abs(d[1]+1) <= 1e-12) {
			return make_rotate(180,Vec3<T>(1,0,0));
		}

		// Householder refl
		Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d R = 2.0 * (d+up) * (d+up).transpose();
		R *= (1.0/denom);
		R -= I;
		for( int i=0; i<3; ++i ){
			for( int j=0; j<3; ++j ){
				xf(i,j) = R(i,j);
			}
		}

		return xf;
	}

	// Makes a view matrix
	// Usage: XForm<float> v = XForm::make_view(eye, viewdir, Vec3f(0,1,0));
	static inline XForm<T> make_view(const Vec3<T> &eye, const Vec3<T> &dir, const Vec3<T> &up)
	{
		Vec3<T> w = dir*-1.f; w.normalize();
		Vec3<T> u = up.cross(w);
		Vec3<T> v = w.cross(u);
		XForm<T> r;
		r.setIdentity();
		for(size_t i=0; i<3; ++i){
			r.data()[4*i] = u[i];
			r.data()[4*i+1] = v[i];
			r.data()[4*i+2] = w[i];
		}
		r.data()[12] = -eye.dot(u);
		r.data()[13] = -eye.dot(v);
		r.data()[14] = -eye.dot(w);
		return r;
	}

	// Makes a view matrix (from a lookat point)
	// Usage: XForm<float> v = XForm::make_lookat(eye, Vec3f(0,0,0), Vec3f(0,1,0));
	static inline XForm<T> make_lookat(
		const Vec3<T> &eye, const Vec3<T> &point, const Vec3<T> &up)
	{
		Vec3<T> dir = point-eye;
		return make_view(eye,dir,up);
	}

	// Makes a perspective matrix
	// Usage: XForm<float> p = XForm::make_persp(45, width/height, 1e-3f, 1e6f);
	static inline XForm<T> make_persp(T fov_deg, T aspect, T near, T far)
	{
		T fov = fov_deg * M_PI / 180.f;
		T cossinf = std::cos(fov/2.f) / std::sin(fov/2.f);
		XForm<T> r;
		r.setIdentity();
		r.data()[0] = cossinf/aspect;
		r.data()[5] = cossinf;
		r.data()[10] = -(near+far)/(far-near);
		r.data()[14] = -(2.f*near*far)/(far-near);
		r.data()[11] = -1.f;
		r.data()[15] = 0.f;
		return r;
	}

	static inline std::string to_string(const XForm<T> &xf)
	{
		std::stringstream ss;
		ss << xf.data()[0];
		for( int i=1; i<16; ++i ){ ss << ' ' << xf.data()[i]; }
		return ss.str();
	}

	static inline XForm<T> from_string(const std::string &s)
	{
		std::stringstream ss; ss << s;
		XForm<T> result;
		// TODO some testing to make sure there are tokens left
		for( int i=0; i<16; ++i ){ ss >> result.data()[i]; }
		return result;
	}

	template<typename U>
	void apply_row(Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic> &V, int row)
	{
		int cols = std::min(dim, (int)V.cols());
		VecD vi = VecD::Zero();
		for (int j=0; j<cols; ++j) { vi[j] = V(row,j); }
		vi = (*this) * vi;
		for (int j=0; j<cols; ++j) { V(row,j) = vi[j]; }
	}

	// Multiply every row (vertex) by the xform
	template<typename U>
	void apply(Eigen::Matrix<U,Eigen::Dynamic,Eigen::Dynamic> &V)
	{
		int nv = V.rows();
		int cols = std::min(dim, (int)V.cols());
		for (int i=0; i<nv; ++i)
		{
			apply_row(V, i);
		}
	}
};

} // ns mcl

#endif
