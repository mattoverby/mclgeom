
// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_EIGENCEREAL_HPP
#define MCL_EIGENCEREAL_HPP 1

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

// From https://stackoverflow.com/questions/22884216/serializing-eigenmatrix-using-cereal-library
// for cereal library (serialization): https://uscilab.github.io/cereal
namespace cereal
{
	template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> static inline
	typename std::enable_if<traits::is_output_serializable<BinaryData<_Scalar>, Archive>::value, void>::type
	save(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> const & m)
	{
		int32_t rows = m.rows();
		int32_t cols = m.cols();
		ar(rows);
		ar(cols);
		if (rows == 0 || cols == 0) { return; }
		ar(binary_data(m.data(), rows * cols * sizeof(_Scalar)));
	}

	template <class Archive, class _Scalar, int _Rows, int _Cols, int _Options, int _MaxRows, int _MaxCols> static inline
	typename std::enable_if<traits::is_input_serializable<BinaryData<_Scalar>, Archive>::value, void>::type
	load(Archive & ar, Eigen::Matrix<_Scalar, _Rows, _Cols, _Options, _MaxRows, _MaxCols> & m)
	{
		int32_t rows = 0;
		int32_t cols = 0;
		ar(rows);
		ar(cols);
		m.resize(rows, cols);
		if (rows == 0 || cols == 0) { return; }
		ar(binary_data(m.data(), static_cast<std::size_t>(rows * cols * sizeof(_Scalar))));
	}

	template <class Archive, class _Scalar, int _Options, typename _StorageIndex> static inline
	typename std::enable_if<traits::is_output_serializable<BinaryData<_Scalar>, Archive>::value, void>::type
	save(Archive & ar, Eigen::SparseMatrix<_Scalar, _Options, _StorageIndex> const & m)
	{
		int32_t nrows = m.rows();
		int32_t ncols = m.cols();
		int32_t nouter = m.outerSize();
		int32_t nnz = m.nonZeros();
		ar(nrows, ncols, nouter, nnz);
		if (nouter == 0 || nnz == 0) { return; }
		std::vector<int32_t> rows, cols;
		std::vector<_Scalar> vals;
		for (int32_t i=0; i<nouter; ++i)
		{
			for(typename Eigen::SparseMatrix<_Scalar,_Options,_StorageIndex>::InnerIterator it(m,i); it; ++it)
			{
				rows.emplace_back(it.row());
				cols.emplace_back(it.col());
				vals.emplace_back(it.value());
			}
		}
		ar(rows, cols, vals);
	}

	template <class Archive, class _Scalar, int _Options, typename _StorageIndex> static inline
	typename std::enable_if<traits::is_input_serializable<BinaryData<_Scalar>, Archive>::value, void>::type
	load(Archive & ar, Eigen::SparseMatrix<_Scalar, _Options, _StorageIndex> & m)
	{
		int32_t nrows = 0;
		int32_t ncols = 0;
		int32_t nouter = 0;
		int32_t nnz = 0;
		ar(nrows, ncols, nouter, nnz);
		m.resize(nrows, ncols);
		if (nouter == 0 || nnz == 0) { return; }
		std::vector<int32_t> rows, cols;
		std::vector<_Scalar> vals;
		ar(rows, cols, vals);
		int32_t n_trip = vals.size();
		std::vector<Eigen::Triplet<_Scalar,_StorageIndex> > triplets;
		triplets.reserve(n_trip);
		for (int32_t i=0; i<n_trip; ++i)
		{
			triplets.emplace_back(rows[i], cols[i], vals[i]);
		}
		m.setFromTriplets(triplets.begin(), triplets.end());
	}
}

#endif
