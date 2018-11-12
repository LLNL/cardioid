#include "SymmetricTensor.hh"
#include "Vector.hh"

SymmetricTensor operator+(const SymmetricTensor& a, const SymmetricTensor& b)
{
   SymmetricTensor sum(a);
   sum.a11 += b.a11;
   sum.a22 += b.a22;
   sum.a33 += b.a33;
   sum.a12 += b.a12;
   sum.a13 += b.a13;
   sum.a23 += b.a23;
   return sum;
}

SymmetricTensor operator/(const SymmetricTensor& m, double d)
{
   SymmetricTensor result(m);
   result.a11 /= d;
   result.a22 /= d;
   result.a33 /= d;
   result.a12 /= d;
   result.a13 /= d;
   result.a23 /= d;
   return result;
}


Vector operator*(const SymmetricTensor& m, const Vector& v)
{
   Vector result;
   result[0] = m.a11*v[0] + m.a12*v[1] + m.a13*v[2];
   result[1] = m.a12*v[0] + m.a22*v[1] + m.a23*v[2];
   result[2] = m.a13*v[0] + m.a23*v[1] + m.a33*v[2];
   
   return result;
}

