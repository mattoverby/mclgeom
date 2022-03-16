// Copyright Matt Overby 2022.
// Distributed under the MIT License.

#ifndef MCL_EIGENMAP_HPP
#define MCL_EIGENMAP_HPP 1

#include <Eigen/Dense>
#include <stdexcept>

namespace mcl
{

template <typename DerivedM, typename DerivedV>
static inline void mat_to_vec(
    const Eigen::PlainObjectBase<DerivedM> &M,
    Eigen::PlainObjectBase<DerivedV> &V)
{
    typedef typename Eigen::PlainObjectBase<DerivedM>::Scalar Scalar;
    int nr = M.rows();
    int nc = M.cols();
    if (V.size() != nr*nc) {
        V.resize(nr*nc);
    }

    if (M.derived().IsRowMajor)
    {
        V = Eigen::Map<DerivedV>((Scalar*)M.data(), M.size());
        return;
    }

    for (int c=0; c<nc; ++c) {
        for (int r=0; r<nr; ++r) {
            V[r*nc+c] = M(r,c);
        }
    }
}

// Matrix must be sized ahead of time
template <typename DerivedM, typename DerivedV>
static inline void vec_to_mat(
    const Eigen::PlainObjectBase<DerivedV> &V,
    Eigen::PlainObjectBase<DerivedM> &M)
{
    typedef typename Eigen::PlainObjectBase<DerivedV>::Scalar Scalar;
    if (M.cols()<1 || M.size() != V.size()) {
        throw std::runtime_error("mcl::vec_to_mat: M not correct size");
    }

    if (M.derived().IsRowMajor)
    {
        M = Eigen::Map<DerivedM>((Scalar*)V.data(), M.rows(), M.cols());
        return;
    } 

    // Column major, need loop.
    int nr = M.rows();
    int nc = M.cols();
    for (int c=0; c<nc; ++c) {
        for (int r=0; r<nr; ++r) {
            M(r,c) = V[r*nc+c];
        }
    }
}

} // ns mcl

#endif