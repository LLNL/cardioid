
#include <simdops/x86_avx512f.hpp>
#include <simdops/x86_avx2.hpp>
#include <simdops/power_qpx.hpp>
#include <simdops/null.hpp>

namespace simdops {

class float64v {
 public:
   inline float64v() {};
   inline float64v(const native_vector_type d) : _d(d) {}
#ifndef SIMDOPS_ARCH_NULL
   explicit inline float64v(const double e) { _d = make_float(e); }
   inline float64v& operator=(const double e) { _d = make_float(e); return *this; }
#endif
   inline float64v& operator=(const native_vector_type d) { _d = d; return *this; }
   
   inline operator native_vector_type() const {  return _d; }
   native_vector_type _d;
};

#ifndef SIMDOPS_ARCH_NULL

inline float64v operator-(const float64v a) { return neg(a); }
//inline float64v operator!(const float64v a) { return b_not(a); }

#define SIMDOPS_BINARY_OP(op,name)                                      \
inline float64v operator op(const float64v a, const float64v b) { return name(a,b); } \
inline float64v operator op(const double a, const float64v b) { return float64v(a) op b; } \
inline float64v operator op(const float64v a, const double b) { return a op float64v(b); }
SIMDOPS_BINARY_OP(+,add)
SIMDOPS_BINARY_OP(*,mul)
SIMDOPS_BINARY_OP(-,sub)
SIMDOPS_BINARY_OP(/,div)
SIMDOPS_BINARY_OP(<,lt)
SIMDOPS_BINARY_OP(>,gt)
SIMDOPS_BINARY_OP(<=,le)
SIMDOPS_BINARY_OP(>=,ge)
SIMDOPS_BINARY_OP(==,eq)
SIMDOPS_BINARY_OP(!=,neq)
SIMDOPS_BINARY_OP(&&,b_and)
SIMDOPS_BINARY_OP(||,b_or)
#undef SIMDOPS_BINARY_OP

#endif

#define SIMDOPS_ASSIGN_OP(op)                                           \
inline float64v& operator op##=(float64v& a, const float64v b) { return (a = a op b); } \
inline float64v& operator op##=(float64v& a, const double b) { return (a = a op b); }
SIMDOPS_ASSIGN_OP(+)
SIMDOPS_ASSIGN_OP(-)
SIMDOPS_ASSIGN_OP(*)
SIMDOPS_ASSIGN_OP(/)
#undef SIMDOPS_ASSIGN_OP

inline double   ternary_if(const double   mask, const double   tt, const double   ff) { return mask ? tt : ff; } 
inline float64v ternary_if(const double   mask, const float64v tt, const float64v ff) { return mask ? tt : ff; } 
inline float64v ternary_if(const double   mask, const double   tt, const float64v ff) { return mask ? float64v(tt) : ff; } 
inline float64v ternary_if(const double   mask, const float64v tt, const double   ff) { return mask ? tt : float64v(ff); }
inline float64v ternary_if(const float64v mask, const float64v tt, const float64v ff) { return b_or(b_and(mask,tt),b_and(b_not(mask),ff)); }
inline float64v ternary_if(const float64v mask, const float64v tt, const double   ff) { return b_or(b_and(mask,tt),b_and(b_not(mask),float64v(ff))); }
inline float64v ternary_if(const float64v mask, const double   tt, const double   ff) { return b_or(b_and(mask,float64v(tt)),b_and(b_not(mask),float64v(ff))); }
inline float64v ternary_if(const float64v mask, const double   tt, const float64v ff) { return b_or(b_and(mask,float64v(tt)),b_and(b_not(mask),ff)); }

inline bool all(const float64v mask) { return !any(b_not(mask)); }
}

#define SIMDOPS_ALIGN(width)

#include <simdops/default_math.hpp>
