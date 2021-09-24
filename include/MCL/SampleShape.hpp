// Copyright Matt Overby 2021.
// Distributed under the MIT License.

#ifndef MCL_SAMPLESHAPE_HPP
#define MCL_SAMPLESHAPE_HPP 1

#include <math.h>

namespace mcl 
{

//
//	Uniform Cone
//
// u1, u2: 0 to 1
static inline void sample_uniform_cone(float u1, float u2, float max_theta, float *vec_3f)
{
  float cos_theta = (1.f - u1) + u1 * cos(max_theta);
  float sin_theta = sqrtf(1.f - cos_theta*cos_theta);
  float phi = u2 * 2.f * M_PI;
  vec_3f[0]=cosf(phi)*sin_theta;
  vec_3f[1]=sinf(phi)*sin_theta;
  vec_3f[2]=cos_theta;
}

//
//	Cosine Hemisphere
//
static inline void sample_cosine_hemisphere(float u1, float u2, float *vec_3f)
{
  float r = sqrt( u1 );
  float theta = 2.f * M_PI * u2;
  vec_3f[0] = r * cosf(theta);
  vec_3f[1] = r * sinf(theta);
  vec_3f[2] = sqrt( fmaxf(0.f, 1.f-u1) );
}

} // end namespace mcl

#endif
