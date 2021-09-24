// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_SIGNEDSVD_HPP
#define MCL_SIGNEDSVD_HPP 1

#include <Eigen/Dense>
#include <Eigen/SVD>

namespace mcl
{

template <typename T, int DIM>
static inline void signed_svd(
	const Eigen::Matrix<T,DIM,DIM> &F,
	Eigen::Matrix<T,DIM,1> &S,
	Eigen::Matrix<T,DIM,DIM> &U,
	Eigen::Matrix<T,DIM,DIM> &V)
{
	using namespace Eigen;
	int dim = DIM == Eigen::Dynamic ? F.rows() : DIM;
	typedef Matrix<T,DIM,DIM> MatX;

	JacobiSVD<MatX> svd(F, ComputeFullU | ComputeFullV);
	S = svd.singularValues();
	U = svd.matrixU();
	V = svd.matrixV();

	MatX J = MatX::Identity(dim,dim);
	J(dim-1,dim-1) = -1.0;

	if (U.determinant() < 0.0)
	{
		U = U * J;
		S[dim-1] *= -1.0;
	}

	if (V.determinant() < 0.0)
	{
		V = (J * V.transpose()).transpose();
		S[dim-1] *= -1.0;
	}

	// Degenerate case
	for(int i=0; i<dim; ++i)
	{
		if (!std::isfinite(S[i]))
		{
			S.setZero();
			U.setIdentity();
			V.setIdentity();
			break;
		}
	}

} // end signed svd

} // end mcl

#endif
