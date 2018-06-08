#pragma once

#if defined(SIMDOPS_ARCH_X86_AVX512F)
#include <immintrin.h>


namespace simdops {
   
typedef __m512d native_vector_type;

#define SIMDOPS_SIZE 8

#define load(x) _mm512_load_pd(x)
#define store(x,y) _mm512_store_pd(x,y)
#define make_float(x) _mm512_set_pd(x,x,x,x,x,x,x,x)
#define splat(x) _mm512_broadcast_f64x4(_mm256_broadcast_sd(&x))
#define add _mm512_add_pd
#define sub _mm512_sub_pd
#define mul _mm512_mul_pd
#define div _mm512_div_pd

#if defined(SIMDOPS_INTEL_VECTOR_LIBM)
inline VReal expm1(VReal x) {
  return _mm512_expm1_pd(x);
}

inline VReal log(VReal x) {
  return _mm512_log_pd(x);
}
#define SIMDOPS_MATH_IS_DEFINED
#endif
   
}

#endif
