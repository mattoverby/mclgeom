// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_COMPUTEROOTS_HPP
#define MCL_COMPUTEROOTS_HPP 1

#include <Eigen/Core>
//#include <cmath>
#include <complex>

namespace mcl
{

// Returns the smallest, positive, real root for quadratic
// a x^2 + b x + c
// Returns a negative value if no roots found
// from https://github.com/mike323zyf/BCQN
static inline double compute_quad_roots(double a, double b, double c, double tol=1e-16)
{
	double t = -1;
	if( std::abs(a) <= tol ){ t = -c / b; }
	else {
		double desc = b*b - 4.0 * a*c;
		if( desc > 0.0 ){
			t = (-b - std::sqrt(desc)) / (2.0 * a);
			if (t < 0.0){
				t = (-b + std::sqrt(desc)) / (2.0 * a);
			}
		}
		else{ // linear
			t = -b / (2.0*a);
		}
	}
	return t;
}

// Returns a negative value if no roots found
// from: https://github.com/libigl/libigl/blob/main/include/igl/flip_avoiding_line_search.cpp
static inline double compute_quad_roots_2(double a, double b, double c)
{
	double t1 = 0, t2 = 0;
	const double polyCoefEps = 1e-16;
	const double max_time = 2;
	if (abs(a) > polyCoefEps)
  {
		double delta_in = pow(b, 2) - 4 * a*c;
		if (delta_in < 0) {
			return -1;
		}
		double delta = sqrt(delta_in);
		t1 = (-b + delta) / (2 * a);
		t2 = (-b - delta) / (2 * a);
	}
	else if (abs(b) > polyCoefEps){
		 t1 = t2 = -c / b;
	} else {
		return -1; 
	}
	if (t1 < 0) t1 = max_time;
	if (t2 < 0) t2 = max_time;

	if (!std::isfinite(t1) || !std::isfinite(t2))
		return -1;

	double tmp_n = std::min(t1, t2);
	t1 = std::max(t1, t2); t2 = tmp_n;
	if (t1 > 0) {
		if (t2 > 0) {
			return t2;
		}
		else {
			return t1;
		}
	}
	else {
		return -1;
	}
}

// Returns the smallest, positive, real root for cubic
// a x^3 + b x^2 + c x + d
// Returns a negative value if no roots found
// from https://github.com/mike323zyf/BCQN
static inline double compute_cubic_roots(double a, double b, double c, double d, double tol=1e-16)
{
	double t = -1;
	if(std::abs(a) <= tol){ t = get_quad_roots(b, c, d, tol); }
	else {
		std::complex<double> i(0, 1);
		std::complex<double> delta0(b*b - 3 * a*c, 0);
		std::complex<double> delta1(2 * b*b*b - 9 * a*b*c + 27 * a*a*d, 0);
		std::complex<double> C = pow((delta1 + sqrt(delta1*delta1 - 4.0 * delta0*delta0*delta0)) / 2.0, 1.0 / 3.0);

		std::complex<double> u2 = (-1.0 + sqrt(3.0)*i) / 2.0;
		std::complex<double> u3 = (-1.0 - sqrt(3.0)*i) / 2.0;

		std::complex<double> t1 = (b + C + delta0 / C) / (-3.0*a);
		std::complex<double> t2 = (b + u2*C + delta0 / (u2*C)) / (-3.0*a);
		std::complex<double> t3 = (b + u3*C + delta0 / (u3*C)) / (-3.0*a);

		if ((std::abs(std::imag(t1))<tol) && (std::real(t1)>0))
			t = std::real(t1);
		if ((std::abs(std::imag(t2))<tol) && (std::real(t2)>0) && ((std::real(t2) < t) || (t < 0)))
			t = std::real(t2);
		if ((std::abs(std::imag(t3))<tol) && (std::real(t3)>0) && ((std::real(t3) < t) || (t < 0)))
			t = std::real(t3);
	}
	return t;
}

// src: https://github.com/libigl/libigl/blob/main/include/igl/flip_avoiding_line_search.cpp
static inline double compute_min_pos_root_2D(const Eigen::MatrixXd& x, const Eigen::MatrixXi& E, Eigen::MatrixXd& p, int f)
{
	int v1 = E(f,0);
	int v2 = E(f,1);
	int v3 = E(f,2);
	const double& U11 = x(v1,0);
	const double& U12 = x(v1,1);
	const double& U21 = x(v2,0);
	const double& U22 = x(v2,1);
	const double& U31 = x(v3,0);
	const double& U32 = x(v3,1);
	const double& V11 = p(v1,0);
	const double& V12 = p(v1,1);
	const double& V21 = p(v2,0);
	const double& V22 = p(v2,1);
	const double& V31 = p(v3,0);
	const double& V32 = p(v3,1);
	double a = V11*V22 - V12*V21 - V11*V32 + V12*V31 + V21*V32 - V22*V31;
	double b = U11*V22 - U12*V21 - U21*V12 + U22*V11 - U11*V32 + U12*V31 + U31*V12 - U32*V11 + U21*V32 - U22*V31 - U31*V22 + U32*V21;
	double c = U11*U22 - U12*U21 - U11*U32 + U12*U31 + U21*U32 - U22*U31;
  double root = get_quad_roots_2(a,b,c);
  if (root < 0) { return std::numeric_limits<float>::max(); }
  return root;
}

} // end mcl

#endif
