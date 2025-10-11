// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_SORT_HPP
#define MCL_SORT_HPP 1

namespace mcl {

template<typename T>
static inline void
sort3(T& a, T& b, T& c)
{
    if (a > b) {
        std::swap(a, b);
    }
    if (b > c) {
        std::swap(b, c);
    }
    if (a > b) {
        std::swap(a, b);
    }
}

} // end ns mcl

#endif
