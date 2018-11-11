#pragma once

#include <cmath>

#if !defined(SIMDOPS_MATH_IS_DEFINED)

namespace simdops {

#define expandMathFunc(name)                            \
inline float64v name(const float64v x)                  \
{                                                       \
   double el[SIMDOPS_FLOAT64V_WIDTH];                   \
   store(el, x);                                        \
   for (int ii=0; ii<SIMDOPS_FLOAT64V_WIDTH; ii++) {    \
      el[ii] = std::name(el[ii]);                       \
   }                                                    \
   float64v whole = load(el);                           \
   return whole;                                        \
}

expandMathFunc(exp)
//expandMathFunc(expm1)
expandMathFunc(log)
expandMathFunc(sqrt)
#undef expandMathFunc

inline float64v pow(const float64v x, double y)
{
   double el[SIMDOPS_FLOAT64V_WIDTH];
   store(el, x);
   for (int ii=0; ii<SIMDOPS_FLOAT64V_WIDTH; ii++) {
      el[ii] = std::pow(el[ii],y);
   }
   float64v whole = load(el);
   return whole;
}

}

#endif
