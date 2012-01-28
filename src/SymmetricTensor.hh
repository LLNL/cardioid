#ifndef SYMMETRIC_TENSOR_HH
#define SYMMETRIC_TENSOR_HH

struct SymmetricTensor
{
   double a11, a12, a13;
   double      a22, a23; // symmetric matrix. Thus a21 = a12
   double           a33; // symmetric matrix. Thus a31 = a13 and a32 = a23
};

#endif

