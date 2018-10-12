#pragma once

#if defined(SIMDOPS_ARCH_POWER_QPX)
#include <builtins.h>


namespace simdops {
   
typedef vector4double native_vector_type;

#define SIMDOPS_FLOAT64V_WIDTH 4

inline native_vector_type load(const double* x) { return vec_ld(0,const_cast<double*>(x)); }
inline void store(double* x, const native_vector_type y) { vec_st(y,0,x); }
inline native_vector_type make_float(const double x) { return vec_splats(x); }
inline native_vector_type splat(const double* x) { return make_float(*x); }
inline native_vector_type add(const native_vector_type a, const native_vector_type b) { return vec_add(a,b); }
inline native_vector_type sub(const native_vector_type a, const native_vector_type b) { return vec_sub(a,b); }
inline native_vector_type mul(const native_vector_type a, const native_vector_type b) { return vec_mul(a,b); }
inline native_vector_type div(const native_vector_type a, const native_vector_type b) { return vec_swdiv_nochk(a,b); }
inline native_vector_type neg(const native_vector_type a) { return vec_neg(a); }

inline native_vector_type b_not(const native_vector_type a) { return vec_not(a); }
inline native_vector_type lt(const native_vector_type a, const native_vector_type b) { return vec_cmplt(b,a); }
inline native_vector_type gt(const native_vector_type a, const native_vector_type b) { return vec_cmpgt(a,b); }
inline native_vector_type le(const native_vector_type a, const native_vector_type b) { return b_not(gt(a,b)); }
inline native_vector_type ge(const native_vector_type a, const native_vector_type b) { return b_not(le(a,b)); }
inline native_vector_type eq(const native_vector_type a, const native_vector_type b) { return vec_cmpeq(a,b); }
inline native_vector_type neq(const native_vector_type a, const native_vector_type b) { return b_not(eq(a,b)); }
inline native_vector_type b_and(const native_vector_type a, const native_vector_type b) { return vec_and(a,b); }
inline native_vector_type b_or(const native_vector_type a, const native_vector_type b) { return vec_or(a,b); }

inline bool any(const native_vector_type a)
{
   for (int ii=0; ii<4; ii++)
   {
      if (vec_extract(a,ii))
      {
         return true;
      }
      return false;
   }
}

}
#endif
