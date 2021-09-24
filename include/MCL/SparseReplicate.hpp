// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_SPARSEREPLICATE_HPP
#define MCL_SPARSEREPLICATE_HPP 1

#include <Eigen/Sparse>

namespace mcl
{

// Copies values along the diagonal of a sparse matrix.
// For making n x m matrices n*dim x m*dim
template<typename T>
static inline void sparse_replicate(
		const Eigen::SparseMatrix<T> &A,
		int dim,
		Eigen::SparseMatrix<T> &Adim)
{
	using namespace Eigen;
	typedef typename SparseMatrix<T>::InnerIterator InnerIter;
	dim = std::max(1,dim);

	Adim.resize(A.rows()*dim, A.cols()*dim);
	Adim.setZero();
	Adim.reserve(A.nonZeros()*dim);

	int cols = A.cols();
	for (int c=0; c<cols; ++c)
	{
		for (int d=0; d<dim; ++d)
		{
			int col = c*dim+d;
			Adim.startVec(col);
			for (InnerIter itA(A,c); itA; ++itA)
			{
				int row = itA.row()*dim+d;
				Adim.insertBack(row, col) = itA.value();
			}
		}
	}

	Adim.finalize();
	Adim.makeCompressed();

} // end sparse replicate

} // ns mcl

#endif
