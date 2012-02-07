#ifndef SYMMETRIC_TENSOR_HH
#define SYMMETRIC_TENSOR_HH

class Vector;

struct SymmetricTensor
{
   double a11, a12, a13;
   double      a22, a23; // symmetric matrix. Thus a21 = a12
   double           a33; // symmetric matrix. Thus a31 = a13 and a32 = a23
};

SymmetricTensor operator+(const SymmetricTensor& a, const SymmetricTensor& b);
Vector operator*(const SymmetricTensor& m, const Vector& v);
SymmetricTensor operator/(const SymmetricTensor& m, double d);
#endif

