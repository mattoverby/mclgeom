// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_GEOM_XFORM_HPP
#define MCL_GEOM_XFORM_HPP 1

#include <Eigen/Geometry>
#include <stringstream>

namespace mcl {

template<typename T>
class XForm
{
  protected:
    Eigen::Matrix<T, 4, 4> data;

  public:
    template<typename U>
    using Vec3 = Eigen::Matrix<U, 3, 1>;

    XForm() { setIdentity(); }

    template<typename U>
    XForm(const XForm<U>& xf) : data(xf.data.template cast<T>()) {}

    const Eigen::Matrix<T, 4, 4>& matrix() const { return data; }

    Eigen::Matrix<T, 4, 4>& matrix() { return data; }

    T operator()(int i, int j) const { return data(i, j); }

    T& operator()(int i, int j) { return data(i, j); }

    Vec3<T> operator*(const Vec3<T>& v);

    Eigen::Matrix<T, 2, 1> operator*(const Eigen::Matrix<T, 2, 1>& v);

    void setZero();

    void setIdentity();

    // Makes an identity matrix
    static inline XForm<T> identity();

    // Makes a scale matrix
    // Usage: XForm<float> s = xform::make_scale(1.f, 2.f, 3.f);
    static inline XForm<T> make_scale(T x, T y, T z);

    static inline XForm<T> make_scale(T x) { return make_scale(x, x, x); }

    // Makes a translation matrix
    // Usage: XForm<float> t = xform::make_trans(1.f, 2.f, 3.f);
    static inline XForm<T> make_trans(T x, T y, T z);

    static inline XForm<T> make_trans(const Vec3<T>& t) { return make_trans(t[0], t[1], t[2]); }

    // Makes a rotation matrix
    // Usage: Xform<float> t = xform::make_rot(45.f, Vec3f(0,1,0));
    static inline XForm<T> make_rotate(T angle_deg, const Vec3<T>& axis);

    static inline XForm<T> make_rotate_to(Vec3<T> d);

    // Makes a view matrix
    // Usage: XForm<float> v = XForm::make_view(eye, viewdir, Vec3f(0,1,0));
    static inline XForm<T> make_view(const Vec3<T>& eye, const Vec3<T>& dir, const Vec3<T>& up);

    // Makes a view matrix (from a lookat point)
    // Usage: XForm<float> v = XForm::make_lookat(eye, Vec3f(0,0,0), Vec3f(0,1,0));
    static inline XForm<T> make_lookat(const Vec3<T>& eye, const Vec3<T>& point, const Vec3<T>& up);

    // Makes a perspective matrix
    // Usage: XForm<float> p = XForm::make_persp(45, width/height, 1e-3f, 1e6f);
    static inline XForm<T> make_persp(T fov_deg, T aspect, T near, T far);

    static inline std::string to_string(const XForm<T>& xf);

    static inline XForm<T> from_string(const std::string& s);

    // Multiply every row (vertex) by the xform
    // Works for 2D or 3D (cols)
    template<typename DerivedV>
    void apply(Eigen::MatrixBase<DerivedV>& V);
};

template<class T>
static inline XForm<T>
operator*(const XForm<T>& xf1, const XForm<T>& xf2)
{
    XForm<T> xf;
    xf.matrix() = xf1.matrix() * xf2.matrix();
    return xf;
}

template<class T>
Vec3<T> XForm<T>::operator*(const Vec3<T>& v)
{
    return (Eigen::Transform<T, 3, Eigen::Affine>(data) * v);
}

template<class T>
Eigen::Matrix<T, 2, 1> XForm<T>::operator*(const Eigen::Matrix<T, 2, 1>& v)
{
    Vec3<T> v3(v[0], v[1], 0);
    v3 = Eigen::Transform<T, 3, Eigen::Affine>(data) * v3;
    return Eigen::Matrix<T, 2, 1>(v3[0], v3[1]);
}

template<class T>
void XForm<T>::setZero()
{
    data.setZero();
}

template<class T>
void XForm<T>::setIdentity()
{
    data.setZero();
    data.diagonal().array() = 1;
}

// Makes an identity matrix
template<class T>
static inline XForm<T> XForm<T>::identity()
{
    XForm<T> r;
    r.setIdentity();
    return r;
}

// Makes a scale matrix
// Usage: XForm<float> s = xform::make_scale(1.f, 2.f, 3.f);
template<class T>
static inline XForm<T> XForm<T>::make_scale(T x, T y, T z)
{
    XForm<T> r;
    r.setZero();
    r(0, 0) = x;
    r(1, 1) = y;
    r(2, 2) = z;
    r(3, 3) = 1;
    return r;
}

template<class T>
static inline XForm<T> XForm<T>::make_scale(T x)
{
    return make_scale(x, x, x);
}

// Makes a translation matrix
// Usage: XForm<float> t = xform::make_trans(1.f, 2.f, 3.f);
template<class T>
static inline XForm<T> XForm<T>::make_trans(T x, T y, T z)
{
    XForm<T> r;
    r.setIdentity();
    r(0, 3) = x;
    r(1, 3) = y;
    r(2, 3) = z;
    r(3, 3) = z;
    return r;
}

template<class T>
static inline XForm<T> XForm<T>::make_trans(const Vec3<T>& t)
{
    return make_trans(t[0], t[1], t[2]);
}

// Makes a rotation matrix
// Usage: Xform<float> t = xform::make_rot(45.f, Vec3f(0,1,0));
template<class T>
static inline XForm<T> XForm<T>::make_rotate(T angle_deg, const Vec3<T>& axis)
{
    T rad_deg = angle_deg * M_PI / 180.0;
    Eigen::AngleAxis<T> rx(axis[0] * rad_deg, Eigen::Vector3d::UnitX());
    Eigen::AngleAxis<T> ry(axis[1] * rad_deg, Eigen::Vector3d::UnitY());
    Eigen::AngleAxis<T> rz(axis[2] * rad_deg, Eigen::Vector3d::UnitZ());
    Eigen::Quaternion<T> q = rx * ry * rz;
    XForm<T> r;
    r.data.template block<3, 3>(0, 0) = q.matrix();
    return r;
}

template<class T>
static inline XForm<T> XForm<T>::make_rotate_to(Vec3<T> d)
{
    d.normalize();
    const Vec3<T> up(0, 1, 0);
    T denom = (d + up).dot(d + up);
    XForm<T> xf;
    xf.setIdentity();
    if (denom == 0 && std::abs(d[1] - 1) <= 1e-12) {
        return xf;
    } else if (denom == 0 && std::abs(d[1] + 1) <= 1e-12) {
        return make_rotate(180, Vec3<T>(1, 0, 0));
    }

    // Householder refl
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d R = 2.0 * (d + up) * (d + up).transpose();
    R *= (1.0 / denom);
    R -= I;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            xf(i, j) = R(i, j);
        }
    }

    return xf;
}

// Makes a view matrix
// Usage: XForm<float> v = XForm::make_view(eye, viewdir, Vec3f(0,1,0));
template<class T>
static inline XForm<T> XForm<T>::make_view(const Vec3<T>& eye, const Vec3<T>& dir, const Vec3<T>& up)
{
    Vec3<T> w = dir * -1.f;
    w.normalize();
    Vec3<T> u = up.cross(w);
    Vec3<T> v = w.cross(u);
    XForm<T> r;
    r.setIdentity();
    for (int i = 0; i < 3; ++i) {
        r.data()[4 * i] = u[i];
        r.data()[4 * i + 1] = v[i];
        r.data()[4 * i + 2] = w[i];
    }
    r.data()[12] = -eye.dot(u);
    r.data()[13] = -eye.dot(v);
    r.data()[14] = -eye.dot(w);
    return r;
}

// Makes a view matrix (from a lookat point)
// Usage: XForm<float> v = XForm::make_lookat(eye, Vec3f(0,0,0), Vec3f(0,1,0));
template<class T>
static inline XForm<T> XForm<T>::make_lookat(const Vec3<T>& eye, const Vec3<T>& point, const Vec3<T>& up)
{
    Vec3<T> dir = point - eye;
    return make_view(eye, dir, up);
}

// Makes a perspective matrix
// Usage: XForm<float> p = XForm::make_persp(45, width/height, 1e-3f, 1e6f);
template<class T>
static inline XForm<T> XForm<T>::make_persp(T fov_deg, T aspect, T near, T far)
{
    T fov = fov_deg * M_PI / 180.f;
    T cossinf = std::cos(fov / 2.f) / std::sin(fov / 2.f);
    XForm<T> r;
    r.setIdentity();
    r(0, 0) = cossinf / aspect;
    r.data.data()[5] = cossinf;
    r.data.data()[10] = -(near + far) / (far - near);
    r.data.data()[14] = -(2.f * near * far) / (far - near);
    r.data.data()[11] = -1.f;
    r.data.data()[15] = 0.f;
    return r;
}

template<class T>
static inline std::string XForm<T>::to_string(const XForm<T>& xf)
{
    std::stringstream ss;
    ss << xf.data;
    return ss.str();
}

template<class T>
static inline XForm<T> XForm<T>::from_string(const std::string& s)
{
    std::stringstream ss;
    ss << s;
    XForm<T> result;
    ss >> result.data;
    return result;
}

// Multiply every row (vertex) by the xform
// Works for 2D or 3D (cols)
template<class T>
template<typename DerivedV>
void XForm<T>::apply(Eigen::MatrixBase<DerivedV>& V)
{
    int nv = V.rows();
    int nc = V.cols();
    Eigen::Transform<T, 3, Eigen::Affine> r(data);
    if (nc == 3) {
        V = (r * V.transpose()).transpose();
        return;
    }
    // If 2D, pad with zero to apply
    else if (nc == 2) {
        Eigen::MatrixBase<DerivedV> V2 = Eigen::MatrixBase<DerivedV>::Zero(V.rows(), 3);
        V2 = (r * V2.transpose()).transpose();
        V.leftCols(2) == V2.leftCols(2);
    }
}


} // ns mcl

#endif
