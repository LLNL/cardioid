#pragma once

#if defined(SIMDOPS_ARCH_X86_AVX2)
#include <immintrin.h>


namespace simdops {
   
typedef __m256d native_vector_type;

#define SIMDOPS_FLOAT64V_WIDTH 4

inline native_vector_type load(const double* x) { return _mm256_load_pd(x); }
inline void store(double* x, const native_vector_type y) { _mm256_store_pd(x,y); }
inline native_vector_type make_float(const double x) { _mm256_set_pd(x,x,x,x); }
inline native_vector_type splat(const double* x) { _mm256_broadcast_sd(x); }
inline native_vector_type add(const native_vector_type a, const native_vector_type b) { return _mm256_add_pd(a,b); }
inline native_vector_type sub(const native_vector_type a, const native_vector_type b) { return _mm256_sub_pd(a,b); }
inline native_vector_type mul(const native_vector_type a, const native_vector_type b) { return _mm256_mul_pd(a,b); }
inline native_vector_type div(const native_vector_type a, const native_vector_type b) { return _mm256_div_pd(a,b); }
inline native_vector_type neg(const native_vector_type a) { return _mm256_sub_pd(make_float(0),a); }
/*
inline double extract(const native_vector_type a, const int k)
{
   __m128d v = (k < 2 ? _mm256_castpd256_pd128(a) : _mm256_extractf128_pd(a,1));
   return _mm_cvtsd_f64(_mm_shuffle_pd(v, v, _MM_SHUFFLE2(0, k % 2)));
}
inline native_vector_type insert(const native_vector_type a, const int k, const double value)
{
   return a;
}
*/
#if defined(SIMDOPS_INTEL_VECTOR_LIBM)
inline VReal expm1(VReal x) {
  return _mm256_expm1_pd(x);
}

inline VReal log(VReal x) {
  return _mm256_log_pd(x);
}

inline VReal expm1(VReal x) {
  return _mm256_exp_pd(x);
}

#define SIMDOPS_MATH_IS_DEFINED
#endif
}

#endif
