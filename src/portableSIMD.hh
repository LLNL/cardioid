#ifndef PORTABLE_SIMD_HH
#define PORTABLE_SIMD_HH

#include <inttypes.h>

// Don't #include this file. #include "ibmInstrinsics.hh" instead.  That
// file contains all of the portability logic so that you don't need any
// #ifdefs to guard the #include.
#ifndef USE_PORTABLE_SIMD
#error "Don't include portableSIMD.hh directly.  Include ibmIntrinsics.hh instead"
#endif

struct vector4double
{
 double v[4]; 
}; 

inline vector4double vec_swdiv_nochk(vector4double op1, vector4double op2)
{
  vector4double target;
  target.v[0]= op1.v[0]/op2.v[0];
  target.v[1]= op1.v[1]/op2.v[1];
  target.v[2]= op1.v[2]/op2.v[2];
  target.v[3]= op1.v[3]/op2.v[3];
  return target;
}

inline vector4double vec_ld(uint32_t shift, float* addr)
{
  vector4double tmp;
  shift /= sizeof(float);
  addr = (float*)(((long long unsigned)addr)  & 0xfffffffffffffff0);
  tmp.v[0] = *(addr + shift);
  tmp.v[1] = *(addr + shift + 1);
  tmp.v[2] = *(addr + shift + 2);
  tmp.v[3] = *(addr + shift + 3);
  return tmp;
}

inline vector4double vec_ld(uint32_t shift, double* addr)
{
  vector4double tmp;
  shift /= sizeof(double);
  addr = (double*)(((long long unsigned)addr)  & 0xffffffffffffffe0);
  tmp.v[0] = *(addr + shift);
  tmp.v[1] = *(addr + shift + 1);
  tmp.v[2] = *(addr + shift + 2);
  tmp.v[3] = *(addr + shift + 3);
  return tmp;
}

inline vector4double vec_splats(double src)
{
  vector4double tmp;
  tmp.v[0] = src;
  tmp.v[1] = src;
  tmp.v[2] = src;
  tmp.v[3] = src;
  return tmp;
}

inline vector4double vec_madd(vector4double op1, vector4double op2, vector4double op3)
{
  vector4double target;
  target.v[0]= op1.v[0]*op2.v[0] + op3.v[0];
  target.v[1]= op1.v[1]*op2.v[1] + op3.v[1];
  target.v[2]= op1.v[2]*op2.v[2] + op3.v[2];
  target.v[3]= op1.v[3]*op2.v[3] + op3.v[3];
  return target;
}

inline vector4double vec_mul(vector4double op1, vector4double op2)
{
  vector4double target;
  target.v[0]= op1.v[0]*op2.v[0];
  target.v[1]= op1.v[1]*op2.v[1];
  target.v[2]= op1.v[2]*op2.v[2];
  target.v[3]= op1.v[3]*op2.v[3];
  return target;
}

inline vector4double vec_add(vector4double op1, vector4double op2)
{
  vector4double target;
  target.v[0]= op1.v[0]+op2.v[0];
  target.v[1]= op1.v[1]+op2.v[1];
  target.v[2]= op1.v[2]+op2.v[2];
  target.v[3]= op1.v[3]+op2.v[3];
  return target;
}

inline vector4double vec_sldw(vector4double op1, vector4double op2, uint32_t shift)
{
  vector4double target;
  target.v[0]= op1.v[shift];
  target.v[1]= ((shift+1) > 3) ? op2.v[shift-3] : op1.v[shift+1];
  target.v[2]= ((shift+2) > 3) ? op2.v[shift-2] : op1.v[shift+2];
  target.v[3]= ((shift+3) > 3) ? op2.v[shift-1] : op1.v[shift+3];
  return target;
}

inline void vec_st(vector4double op1, uint32_t offset, double* addr)
{
  offset /= sizeof(double);
  addr  = (double *)(((long long unsigned)addr)  & 0xffffffffffffffe0);
  *(addr + offset)     = op1.v[0];
  *(addr + offset + 1) = op1.v[1];
  *(addr + offset + 2) = op1.v[2];
  *(addr + offset + 3) = op1.v[3];
}

#endif
