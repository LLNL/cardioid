#pragma once

#if defined(SIMDOPS_ARCH_X86_AVX512F)
#include <immintrin.h>


namespace simdops {
   
typedef __m512d native_vector_type;

#define SIMDOPS_FLOAT64V_WIDTH 8

inline native_vector_type load(const double* x) { return _mm512_load_pd(x); }
inline void store(double* x, const native_vector_type y) { _mm512_store_pd(x,y); }
inline native_vector_type make_float(const double x) { return _mm512_set_pd(x,x,x,x,x,x,x,x); }
inline native_vector_type splat(const double* x) { return _mm512_broadcast_f64x4(_mm256_broadcast_sd(x)); }
inline native_vector_type add(const native_vector_type a, const native_vector_type b) { return _mm512_add_pd(a,b); }
inline native_vector_type sub(const native_vector_type a, const native_vector_type b) { return _mm512_sub_pd(a,b); }
inline native_vector_type mul(const native_vector_type a, const native_vector_type b) { return _mm512_mul_pd(a,b); }
inline native_vector_type div(const native_vector_type a, const native_vector_type b) { return _mm512_div_pd(a,b); }
inline native_vector_type neg(const native_vector_type a) { return _mm512_sub_pd(make_float(0),a); }

#if defined(SIMDOPS_INTEL_VECTOR_LIBM)
inline native_vector_type expm1(native_vector_type x) {
  return _mm512_expm1_pd(x);
}

inline native_vector_type log(native_vector_type x) {
  return _mm512_log_pd(x);
}
   
inline native_vector_type exp(native_vector_type x) {
  return _mm512_exp_pd(x);
}
#define SIMDOPS_MATH_IS_DEFINED
#endif
   
}

#endif
