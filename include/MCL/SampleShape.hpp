// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_SAMPLESHAPE_HPP
#define MCL_SAMPLESHAPE_HPP 1

#include <math.h>

namespace mcl {

//
//	Uniform Cone
//
// u1, u2: 0 to 1
template<typename T>
static inline void
sample_uniform_cone(T u1, T u2, T max_theta, T* vec_3)
{
    T cos_theta = (T(1) - u1) + u1 * std::cos(max_theta);
    T sin_theta = std::sqrt(1 - cos_theta * cos_theta);
    T phi = u2 * 2 * M_PI;
    vec_3[0] = std::cos(phi) * sin_theta;
    vec_3[1] = std::sin(phi) * sin_theta;
    vec_3[2] = cos_theta;
}

//
//	Cosine Hemisphere
//
template<typename T>
static inline void
sample_cosine_hemisphere(T u1, T u2, T* vec_3)
{
    T r = std::sqrt(u1);
    T theta = 2 * M_PI * u2;
    vec_3[0] = r * std::cos(theta);
    vec_3[1] = r * std::sin(theta);
    vec_3[2] = std::sqrt(std::max(T(0), 1 - u1));
}

} // end namespace mcl

#endif
