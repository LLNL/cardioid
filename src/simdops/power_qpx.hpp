#pragma once

#if defined(SIMDOPS_ARCH_POWER_QPX)
#include <builtins.h>


namespace simdops {
   
typedef vector4double native_vector_type;

#define SIMDOPS_FLOAT64V_WIDTH 4

inline native_vector_type load(const double* x) { return vec_ld(x); }
inline void store(double* x, const native_vector_type y) { vec_store(y,0,x); }
inline native_vector_type make_float(const double x) { return vec_spats(x); }
inline native_vector_type splat(const double* x) { return make_float(*x); }
inline native_vector_type add(const native_vector_type a, const native_vector_type b) { return vec_add(a,b); }
inline native_vector_type sub(const native_vector_type a, const native_vector_type b) { return vec_sub(a,b); }
inline native_vector_type mul(const native_vector_type a, const native_vector_type b) { return vec_mul(a,b); }
inline native_vector_type div(const native_vector_type a, const native_vector_type b) { return vec_swdiv_nochk(a,b); }
inline native_vector_type neg(const native_vector_type a) { return vec_neg(a); }

inline native_vector_type lt(const native_vector_type a, const native_vector_type b) { return vec_cmpgt(b,a); }
inline native_vector_type gt(const native_vector_type a, const native_vector_type b) { return vec_cmpgt(a,b); }
inline native_vector_type le(const native_vector_type a, const native_vector_type b) { return vec_cmpge(b,a); }
inline native_vector_type ge(const native_vector_type a, const native_vector_type b) { return vec_cmpge(a,b); }
inline native_vector_type eq(const native_vector_type a, const native_vector_type b) { return vec_cmpeq(a,b); }
inline native_vector_type neq(const native_vector_type a, const native_vector_type b) { return b_not(eq(a,b)); }
inline native_vector_type b_and(const native_vector_type a, const native_vector_type b) { return vec_and(a,b); }
inline native_vector_type b_or(const native_vector_type a, const native_vector_type b) { return vec_or(a,b); }
inline native_vector_type b_not(const native_vector_type a) { return vec_not(a); }

}

#endif
