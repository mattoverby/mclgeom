
// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_EIGENCEREAL_HPP
#define MCL_EIGENCEREAL_HPP 1

#include <Eigen/Core>
#include <cereal/archives/binary.hpp>

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
		int32_t rows;
		int32_t cols;
		ar(rows);
		ar(cols);
		m.resize(rows, cols);
		if (rows == 0 || cols == 0) { return; }
		ar(binary_data(m.data(), static_cast<std::size_t>(rows * cols * sizeof(_Scalar))));
	}
}

#endif
