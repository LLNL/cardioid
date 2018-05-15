
#include <simdops/x86_avx512f.hpp>
#include <simdops/x86_avx2.hpp>
#include <simdops/power_qpx.hpp>
#include <simdops/null.hpp>

namespace simdops {

class float64v {
 public:
   inline float64v() {};
   inline float64v(const native_vector_type d) : _d(d) {}
   explicit inline float64v(const double e) { _d = make_float(e); }
   inline float64v& operator=(const native_vector_type d) { _d = d; return *this; }
   inline operator native_vector_type() const {  return _d; }
   native_vector_type _d;
};
   
inline float64v operator+(const float64v a, const float64v b)
{
   return add(a,b);
}
inline float64v operator*(const float64v a, const float64v b)
{
   return mul(a,b);
}
inline float64v operator-(const float64v a, const float64v b)
{
   return sub(a,b);
}
inline float64v operator/(const float64v a, const float64v b)
{
   return div(a,b);
}
inline float64v operator-(const float64v a)
{
   return neg(a);
}

inline float64v operator+(const double a, const float64v b)
{
  return float64v(a)+b;
}
inline float64v operator+(const float64v a, const double b)
{
  return a+float64v(b);
}

inline float64v operator*(const double a, const float64v b)
{
  return float64v(a)*b;
}
inline float64v operator*(const float64v a, const double b)
{
  return a*float64v(b);
}

inline float64v operator-(const double a, const float64v b)
{
  return float64v(a)-b;
}
inline float64v operator-(const float64v a, const double b)
{
  return a-float64v(b);
}

inline float64v operator/(const double a, const float64v b)
{
  return float64v(a)/b;
}
inline float64v operator/(const float64v a, const double b)
{
  return a/float64v(b);
}

inline float64v& operator+=(float64v& a, const float64v b)
{
   return (a = a+b);
}

inline float64v& operator-=(float64v& a, const float64v b)
{
   return (a = a-b);
}

inline float64v& operator*=(float64v& a, const float64v b)
{
   return (a = a*b);
}

inline float64v& operator/=(float64v& a, const float64v b)
{
   return (a = a/b);
}

}

#define SIMDOPS_ALIGN(width)

#include <simdops/default_math.hpp>
