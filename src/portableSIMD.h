#ifndef PORTABLE_SIMD_HH
#define PORTABLE_SIMD_HH

// use the portable replacement funtions by default.
#define USE_PORTABLE_SIMD


// BGQ has the real SIMD functions.  Don't use replacements
#ifdef BGQ
#undef USE_PORTABLE_SIMD
#endif

#ifdef USE_PORTABLE_SIMD


#ifdef __cplusplus
extern "C" {
#endif

#include <inttypes.h>

/*
Don't %include this file. %include "ibmInstrinsics.hh" instead.  That
file contains all of the portability logic so that you don't need any
%ifdefs to guard the %include.
*/
#ifndef USE_PORTABLE_SIMD
#error "Don't include portableSIMD.hh directly.  Include ibmIntrinsics.hh instead"
#endif
#ifndef unint32_t
#define uint32_t unsigned 
#endif

typedef struct vector4double_
{
 double v[4]; 
} vector4double; 

#define vec_ld(shift,addr)  vload((shift),(void*)(addr),sizeof((addr))) 

vector4double vload(uint32_t shift, void* addr, int size)
{
  vector4double tmp;
  shift /= size;
  if (size==4) 
  {
  tmp.v[0] = *(((float*)addr) + shift);
  tmp.v[1] = *(((float*)addr) + shift + 1);
  tmp.v[2] = *(((float*)addr) + shift + 2);
  tmp.v[3] = *(((float*)addr) + shift + 3);
  return tmp;
  }
  if (size==8) 
  {
  tmp.v[0] = *(((double*)addr) + shift);
  tmp.v[1] = *(((double*)addr) + shift + 1);
  tmp.v[2] = *(((double*)addr) + shift + 2);
  tmp.v[3] = *(((double*)addr) + shift + 3);
  return tmp;
  }
}
   
/*
vector4double vec_ld(uint32_t shift, double* addr)
{
  vector4double tmp;
  shift /= sizeof(double);
  tmp.v[0] = *(addr + shift);
  tmp.v[1] = *(addr + shift + 1);
  tmp.v[2] = *(addr + shift + 2);
  tmp.v[3] = *(addr + shift + 3);
  return tmp;
}
*/
vector4double vec_lds(uint32_t shift, double* addr)
{
  vector4double tmp;
  shift /= sizeof(double);
  tmp.v[0] = *(addr + shift);
  tmp.v[1] = *(addr + shift);
  tmp.v[2] = *(addr + shift);
  tmp.v[3] = *(addr + shift);
  return tmp;
}

vector4double vec_splats(double src)
{
  vector4double tmp;
  tmp.v[0] = src;
  tmp.v[1] = src;
  tmp.v[2] = src;
  tmp.v[3] = src;
  return tmp;
}

vector4double vec_re(vector4double op)
{
  double eps = 1/4096; 
  vector4double target;
  target.v[0]= (1.0+eps)/(op.v[0]);
  target.v[1]= (1.0-eps)/(op.v[1]);
  target.v[2]= (1.0+eps)/(op.v[2]);
  target.v[3]= (1.0-eps)/(op.v[3]);
  return target; 
}
vector4double vec_swdiv_nochk(vector4double op1, vector4double op2)
{
  vector4double target;
  target.v[0]= op1.v[0]/op2.v[0];
  target.v[1]= op1.v[1]/op2.v[1];
  target.v[2]= op1.v[2]/op2.v[2];
  target.v[3]= op1.v[3]/op2.v[3];
  return target;
}

vector4double vec_madd(vector4double op1, vector4double op2, vector4double op3)
{
  vector4double target;
  target.v[0]= op1.v[0]*op2.v[0] + op3.v[0];
  target.v[1]= op1.v[1]*op2.v[1] + op3.v[1];
  target.v[2]= op1.v[2]*op2.v[2] + op3.v[2];
  target.v[3]= op1.v[3]*op2.v[3] + op3.v[3];
  return target;
}
vector4double vec_nmsub(vector4double op1, vector4double op2, vector4double op3)
{
  vector4double target;
  target.v[0]= -op1.v[0]*op2.v[0] + op3.v[0];
  target.v[1]= -op1.v[1]*op2.v[1] + op3.v[1];
  target.v[2]= -op1.v[2]*op2.v[2] + op3.v[2];
  target.v[3]= -op1.v[3]*op2.v[3] + op3.v[3];
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

vector4double vec_sub(vector4double op1, vector4double op2)
{
  vector4double target;
  target.v[0]= op1.v[0]-op2.v[0];
  target.v[1]= op1.v[1]-op2.v[1];
  target.v[2]= op1.v[2]-op2.v[2];
  target.v[3]= op1.v[3]-op2.v[3];
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

void vec_st(vector4double op1, uint32_t offset, double* addr)
{
  offset /= sizeof(double);
  *(addr + offset)     = op1.v[0];
  *(addr + offset + 1) = op1.v[1];
  *(addr + offset + 2) = op1.v[2];
  *(addr + offset + 3) = op1.v[3];
}

#endif
#endif
