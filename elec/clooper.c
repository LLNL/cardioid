#include "clooper.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

/*
   This file contains four different routines:
     integrateLoop_normal(),
     integrateLoop_vector(),
     integrateLoop_compare(), and
     integrateLoop2()

   All routines have the same calling sequence. The main
   calling code (in simulationLoop.cc) expects a routine
   named integrateLoop(). In order to make it easy to switch
   between the different routines, there is a macro below this
   comment which redefines one of the above routine names to
   integrateLoop. clooper.h declares one function, namely
   integrateLoop().

   Description of rouines:
     integrateLoop_normal() is the original code.

     integrateLoop_vector() is a vectorized routine (for BG/Q)
 
     integrateLoop_compare(), calls both integrateLoop_normal() and
     integrateLoop_vector(), with separate output arrays, and compares
     all relevant output entries for equality. It prints a message and
     halts the program if there is any difference. It can be used as a
     debugging/validation routine to check correctness of the vectorized
     routine.

     integrateLoop2(), seems like a remnant of old days,
     it is basically like integrateLoop(), but unrolled
     by doing two cells per loop iteration.

     Except for integrateLoop2() and the body of integrateLoop_normal(),
     blame T. Oppelstrup for any problems.
*/

#ifdef BGQ
#define integrateLoop_vector integrateLoop
#define integrateLoop_vector_nostim integrateLoop_nostim
#else
#define integrateLoop_normal integrateLoop
#define integrateLoop_normal_nostim integrateLoop_nostim
#endif
/* Other options for the "integrateLoop" macro...
//#define integrateLoop_normal integrateLoop
//#define integrateLoop_compare integrateLoop
*/

void integrateLoop_normal(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale)
{
   for (unsigned ii=begin; ii<end; ++ii)
   {
      int ib = blockOffset[ii];

      double dVm = dVmR[ii] + stim[ii] + diffusionScale * dVmBlock[ib];
      Vm[ii] += dt*dVm;
      VmBlock[ib] = Vm[ii];
   }
   return;
}
void integrateLoop_normal_nostim(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale)
{
   for (unsigned ii=begin; ii<end; ++ii)
   {
      int ib = blockOffset[ii];

      double dVm = dVmR[ii] + /*stim[ii] +*/ diffusionScale * dVmBlock[ib];
      Vm[ii] += dt*dVm;
      VmBlock[ib] = Vm[ii];
   }
   return;
}

void integrateLoopLag(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale)
{
   double dtScale = dt*diffusionScale; 
   for (unsigned ii=begin; ii<end; ++ii)
   {
      int ib = blockOffset[ii];

      VmBlock[ib] = Vm[ii];
      Vm[ii] += dt*dVmR[ii] + dt*stim[ii] + dtScale * dVmBlock[ib];
   }
   return;
}

#ifdef BGQ
void integrateLoop_vector(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale)
{
  unsigned int n = end - begin;
  unsigned ii;

  dVmR += begin - 1;
  stim += begin - 1;
  blockOffset += begin-1;
  Vm += begin - 1;
  
  for(ii = 0; ii<n; ii++) {
    if( (((unsigned long long int) Vm)/sizeof(double)) % 4 == 3) break;
    int ib = *(++blockOffset);

    double dVm = *(++dVmR) + *(++stim) + diffusionScale * dVmBlock[ib];
    double dv = dt*dVm + *(++Vm);
    *Vm = dv;
    VmBlock[ib] = dv;
  }

  n = n - ii;
  if(n > 0) {
    vector4double diffScale_val = {diffusionScale,diffusionScale,diffusionScale,diffusionScale};
    vector4double dt_val = {dt,dt,dt,dt};

#define makeshift(name) \
    vector4double name##_shift = vec_lvsl(0,name-3); \
    vector4double *name##_ptr = (vector4double *) (name-3); \
    vector4double name##_val,name##_v1,name##_v2 = *(++name##_ptr)
#define shiftload(name) \
    name##_v1 = name##_v2; \
    name##_v2 = *(++name##_ptr); \
    name##_val = vec_perm(name##_v1,name##_v2,name##_shift)

    const int n4 = n>>2;
    makeshift(dVmR);
    makeshift(stim);
    vector4double *Vm_ptr = (vector4double *) (Vm-3);
    for(ii = 0; ii<n4; ii++) {
      vector4double Vm_val = *(++Vm_ptr),x;
      shiftload(dVmR);
      shiftload(stim);

      int ib1 = *(++blockOffset) * sizeof(double);
      int ib2 = *(++blockOffset) * sizeof(double);
      int ib3 = *(++blockOffset) * sizeof(double);
      int ib4 = *(++blockOffset) * sizeof(double);

      {
	vector4double y1,y2,y3;
	//x[0] = dVmBlock[ib];
	asm("lfdx %0,%1,%2" : "=v"(x) : "b"(ib1) , "b"(dVmBlock));
	asm("lfdx %0,%1,%2" : "=v"(y1) : "b"(ib2) , "b"(dVmBlock));
	x = vec_sldw(x,y1,1);
	asm("lfdx %0,%1,%2" : "=v"(y2) : "b"(ib3) , "b"(dVmBlock));
	x = vec_sldw(x,y2,1);
	asm("lfdx %0,%1,%2" : "=v"(y3) : "b"(ib4) , "b"(dVmBlock));
	x = vec_sldw(x,y3,1);
      }

      vector4double dVm = vec_madd(diffScale_val , x, stim_val) ;
      dVm = vec_add(dVm, dVmR_val);
      vector4double dv = vec_madd(  dt_val,dVm , Vm_val);
      *Vm_ptr = dv;
      asm("stfdx %0,%1,%2" : : "v"(dv ) , "b"(ib1) , "b"(VmBlock));
      vector4double dv1 = vec_sldw(dv,dv,1);
      asm("stfdx %0,%1,%2" : : "v"(dv1) , "b"(ib2) , "b"(VmBlock));
      vector4double dv2 = vec_sldw(dv,dv,2);
      asm("stfdx %0,%1,%2" : : "v"(dv2) , "b"(ib3) , "b"(VmBlock));
      vector4double dv3 = vec_sldw(dv,dv,3);
      asm("stfdx %0,%1,%2" : : "v"(dv3) , "b"(ib4) , "b"(VmBlock));
    }

    ii = ii << 2;
    if(ii < n) {
      dVmR += ii;
      stim += ii;
      Vm += ii;
      for(; ii<n; ii++) {
	int ib = *(++blockOffset);
	double dVm = *(++dVmR) + *(++stim) + diffusionScale * dVmBlock[ib];
	double dv = dt*dVm + *(++Vm);
	*Vm = dv;
	VmBlock[ib] = dv;
      }
    }

#undef shiftload
#undef makeshift

  }

}

void integrateLoop_vector_nostim(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale)
{
  unsigned int n = end - begin;
  unsigned ii;

  dVmR += begin - 1;
  blockOffset += begin-1;
  Vm += begin - 1;
  
  for(ii = 0; ii<n; ii++) {
    if( (((unsigned long long int) Vm)/sizeof(double)) % 4 == 3) break;
    int ib = *(++blockOffset);
    double dVm = *(++dVmR) + diffusionScale * dVmBlock[ib];
    double dv = dt*dVm + *(++Vm);
    *Vm = dv;
    VmBlock[ib] = dv;
  }

  n = n - ii;
  if(n > 0) {
    vector4double diffScale_val = {diffusionScale,diffusionScale,diffusionScale,diffusionScale};
    vector4double dt_val = {dt,dt,dt,dt};

#define makeshift(name) \
    vector4double name##_shift = vec_lvsl(0,name-3); \
    vector4double *name##_ptr = (vector4double *) (name-3); \
    vector4double name##_val,name##_v1,name##_v2 = *(++name##_ptr)
#define shiftload(name) \
    name##_v1 = name##_v2; \
    name##_v2 = *(++name##_ptr); \
    name##_val = vec_perm(name##_v1,name##_v2,name##_shift)

    const int n4 = n>>2;
    makeshift(dVmR);
    vector4double *Vm_ptr = (vector4double *) (Vm-3);
    for(ii = 0; ii<n4; ii++) {
      vector4double Vm_val = *(++Vm_ptr),x;
      shiftload(dVmR);

      int ib1 = *(++blockOffset) * sizeof(double);
      int ib2 = *(++blockOffset) * sizeof(double);
      int ib3 = *(++blockOffset) * sizeof(double);
      int ib4 = *(++blockOffset) * sizeof(double);

      {
	vector4double y1,y2,y3;
	//x[0] = dVmBlock[ib];
	asm("lfdx %0,%1,%2" : "=v"(x) : "b"(ib1) , "b"(dVmBlock));
	asm("lfdx %0,%1,%2" : "=v"(y1) : "b"(ib2) , "b"(dVmBlock));
	x = vec_sldw(x,y1,1);
	asm("lfdx %0,%1,%2" : "=v"(y2) : "b"(ib3) , "b"(dVmBlock));
	x = vec_sldw(x,y2,1);
	asm("lfdx %0,%1,%2" : "=v"(y3) : "b"(ib4) , "b"(dVmBlock));
	x = vec_sldw(x,y3,1);
      }

      vector4double dVm = vec_madd(diffScale_val , x, dVmR_val) ;
      vector4double dv = vec_madd(  dt_val,dVm , Vm_val);
      *Vm_ptr = dv;
      asm("stfdx %0,%1,%2" : : "v"(dv ) , "b"(ib1) , "b"(VmBlock));
      vector4double dv1 = vec_sldw(dv,dv,1);
      asm("stfdx %0,%1,%2" : : "v"(dv1) , "b"(ib2) , "b"(VmBlock));
      vector4double dv2 = vec_sldw(dv,dv,2);
      asm("stfdx %0,%1,%2" : : "v"(dv2) , "b"(ib3) , "b"(VmBlock));
      vector4double dv3 = vec_sldw(dv,dv,3);
      asm("stfdx %0,%1,%2" : : "v"(dv3) , "b"(ib4) , "b"(VmBlock));
    }

    ii = ii << 2;
    if(ii < n) {
      dVmR += ii;
      stim += ii;
      Vm += ii;
      for(; ii<n; ii++) {
	int ib = *(++blockOffset);
	double dVm = *(++dVmR) + diffusionScale * dVmBlock[ib];
	double dv = dt*dVm + *(++Vm);
	*Vm = dv;
	VmBlock[ib] = dv;
      }
    }

#undef shiftload
#undef makeshift

  }

}
#endif

#ifdef BGQ
void integrateLoop_compare(const int begin, const int end,
                   const double dt,
                   double* dVmR,
                   double* stim,
                   unsigned* blockOffset,
                   double* dVmBlock,
                   double* VmBlock,
                   double* Vm,
                   double diffusionScale)
{
  double *Vm2_base = (double *) malloc(sizeof(double) * (end-begin)),*Vm2 = Vm2_base - begin;
  double *VmBlock2,*VmBlock2_base;
  
  {
    unsigned int vmb0 = 0,vmb1 = 0,i;
    for(i = begin; i<end; i++) {
      unsigned int idx = blockOffset[i];
      if(i == 0 || idx < vmb0) vmb0 = idx;
      if(i == 0 || idx > vmb1) vmb1 = idx;
      Vm2[i] = Vm[i];
    }
    VmBlock2_base = (double *) malloc(sizeof(double) * (vmb1 - vmb0 + 1));
    VmBlock2 = VmBlock2_base - vmb0;
  }


  integrateLoop_normal(begin,end,dt,dVmR,stim,blockOffset,dVmBlock,
		       VmBlock,Vm,diffusionScale);

  integrateLoop_vector(begin,end,dt,dVmR,stim,blockOffset,dVmBlock,
		       VmBlock2,Vm2,diffusionScale);


  {
    int i;
    for(i = begin; i<end; i++) {
      unsigned int ib = blockOffset[i];
      if(Vm[i] != Vm2[i]) {
	printf("%s:%d: Vm error at index %d (%d,%d) :: Vm = %15.5e , Vm2 = %15.5e\n",
	       __FILE__,__LINE__,i,begin,end,Vm[i],Vm2[i]);
	exit(1);
      }
      if(VmBlock[ib] != VmBlock2[ib]) {
	printf("%s:%d: VmBlock error at index %d,%d (%d,%d) :: Vm = %15.5e , Vm2 = %15.5e\n",
	       __FILE__,__LINE__,i,ib,begin,end,VmBlock[ib],VmBlock2[ib]);
	exit(1);
      }
    }
  }

  free(VmBlock2_base);
  free(Vm2_base);
}
#endif

/* This routine was in the old clooper.c, and seems like an attempt to
   pipeline/iunroll the code. It was not used before, and is not used
   now.
*/
void integrateLoop2(const int size, const double dt, double* dVmR, double* dVmD, double* Vm)
{
   const int halfSize = size/2;
   for (unsigned ii=0; ii<halfSize; ++ii)
   {
      double dVm1 = dVmR[ii] + dVmD[ii];
      double dVm2 = dVmR[ii+halfSize] + dVmD[ii+halfSize];
      Vm[ii] += dt*dVm1;
      Vm[ii+halfSize] += dt*dVm2;
   }
   if (size%2 != 0)
   {
      const int jj = size-1;
      double dVm1 = dVmR[jj] + dVmD[jj];
      Vm[jj] += dt*dVm1;
   }
   return;
}
