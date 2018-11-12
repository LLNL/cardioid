#pragma once

#if !defined(SIMDOPS_FLOAT64V_WIDTH)

namespace simdops {
typedef double native_vector_type;

#define SIMDOPS_FLOAT64V_WIDTH 1
#define SIMDOPS_SIZE 1
#define SIMDOPS_WIDTH 1

inline native_vector_type load(const double* x) { return *x; }
inline void store(double* x, const native_vector_type y) { *x = y; }
inline native_vector_type make_float(const double x) { return x; }
inline native_vector_type splat(const double* x) { return *x; }
inline native_vector_type add(const native_vector_type a, const native_vector_type b) { return a+b; }
inline native_vector_type sub(const native_vector_type a, const native_vector_type b) { return a-b; }
inline native_vector_type mul(const native_vector_type a, const native_vector_type b) { return a*b; }
inline native_vector_type div(const native_vector_type a, const native_vector_type b) { return a/b; }
inline native_vector_type neg(const native_vector_type a) { return -a; }

inline native_vector_type lt(const native_vector_type a, const native_vector_type b) { return a<b; }
inline native_vector_type gt(const native_vector_type a, const native_vector_type b) { return a>b; }
inline native_vector_type le(const native_vector_type a, const native_vector_type b) { return a<=b; }
inline native_vector_type ge(const native_vector_type a, const native_vector_type b) { return a>=b; }
inline native_vector_type eq(const native_vector_type a, const native_vector_type b) { return a==b; }
inline native_vector_type neq(const native_vector_type a, const native_vector_type b) { return a!=b; }
inline native_vector_type b_and(const native_vector_type a, const native_vector_type b) { return a && b; }
inline native_vector_type b_or(const native_vector_type a, const native_vector_type b) { return a || b; }
inline native_vector_type b_not(const native_vector_type a) { return !a; }
   
inline bool any(const native_vector_type a) { return a; }

}

   
#define SIMDOPS_MATH_IS_DEFINED
#define SIMDOPS_COMPILER_HAS_DEFAULT_OPS
#define SIMDOPS_ARCH_NULL
#endif
