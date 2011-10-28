#include "Vector.hh"

#include <cmath>
#include <iostream>

using std::sqrt;
using std::ostream;

Vector::Vector()
{
   _r[0]=0.0;
   _r[1]=0.0;
   _r[2]=0.0;
}

Vector::Vector(double x, double y, double z)
{
   _r[0]=x;
   _r[1]=y;
   _r[2]=z;
}
   
void
Vector::assign(double x, double y, double z)
{
   _r[0] = x;
   _r[1] = y;
   _r[2] = z;
}


double
Vector::length() const
{
   return sqrt(_r[0]*_r[0] + _r[1]*_r[1] + _r[2]*_r[2]);
}

void
Vector::normalize()
{
   double len = length();
   if ( len > 0 )
      *this/=len; 
}
   
Vector&
Vector::operator+=(const Vector& a)
{
   this->_r[0] += a._r[0];
   this->_r[1] += a._r[1];
   this->_r[2] += a._r[2];
   return *this;
}
Vector&
Vector::operator-=(const Vector& a)
{
   this->_r[0] -= a._r[0];
   this->_r[1] -= a._r[1];
   this->_r[2] -= a._r[2];
   return *this;
}
Vector&
Vector::operator*=(const double& a)
{
   this->_r[0] *= a;
   this->_r[1] *= a;
   this->_r[2] *= a;
   return *this;
}
Vector&
Vector::operator/=(const double& a)
{
   this->_r[0] /= a;
   this->_r[1] /= a;
   this->_r[2] /= a;
   return *this;
}

Vector
Vector::operator-() const
{
   return (*this) * -1.0;
}


bool
Vector::operator<(const Vector& b) const
{
   return
      _r[0]<b._r[0] ||
      (_r[0]==b._r[0] &&
       (_r[1]<b._r[1] ||
	( _r[1] == b._r[1] && _r[2] < b._r[2] )));
}

Vector
operator+(const Vector& a, const Vector& b)
{
   Vector c(a);
   return c+=b;
}
Vector
operator-(const Vector& a, const Vector& b)
{
   Vector c(a);
   return c-=b;
}
Vector
operator*(const Vector& a, const double& b)
{
   Vector c(a);
   return c*=b;
}
Vector
operator*(const double& a, const Vector& b)
{
   Vector c(b);
   return c*=a;
}
Vector
operator/(const Vector& a, const double& b)
{
   Vector c(a);
   return c/=b;
}
Vector
operator/(const double& a, const Vector& b)
{
   Vector c(b);
   return c/=a;
}
double
dot(const Vector& a, const Vector& b)
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
double
length(const Vector& a)
{
   return a.length();
}

ostream&
operator<<(ostream& os, const Vector& r)
{
   os << "(" << r[0] << ", " << r[1] << ", " << r[2] << ")";
   return os;
}
