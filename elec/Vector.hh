#ifndef VECTOR_HH
#define VECTOR_HH

#include <iosfwd>
#include <cmath>
class Vector
{
 public:

   Vector();
   Vector(double x, double y, double z);
   
   double& operator[](unsigned ii)             {return _r[ii];}
   const double& operator[](unsigned ii) const {return _r[ii];}

   double&       x()       {return _r[0];}
   const double& x() const {return _r[0];}
   double&       y()       {return _r[1];}
   const double& y() const {return _r[1];}
   double&       z()       {return _r[2];}
   const double& z() const {return _r[2];}

   void assign(double x, double y, double z);
   
   double length() const;
   void normalize();
   
   Vector& operator+=(const Vector& a);
   Vector& operator-=(const Vector& a);
   Vector& operator*=(const double& a);
   Vector& operator/=(const double& a);
   Vector operator-() const;
   bool operator<(const Vector& b) const;

 private:
   double _r[3];
};

Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(const Vector& a, const double& b);
Vector operator*(const double& a, const Vector& b);
Vector operator/(const Vector& a, const double& b);
Vector operator/(const double& a, const Vector& b);
double dot(const Vector& a, const Vector& b);
double diffSq(const Vector& a, const Vector& b);
double length(const Vector& a);

// input , output functions
std::ostream& operator<<(std::ostream&, const Vector&);

// Inline Implementations:

inline
Vector::Vector()
{
   _r[0]=0.0;
   _r[1]=0.0;
   _r[2]=0.0;
}

inline
Vector::Vector(double x, double y, double z)
{
   _r[0]=x;
   _r[1]=y;
   _r[2]=z;
}

inline void
Vector::assign(double x, double y, double z)
{
   _r[0] = x;
   _r[1] = y;
   _r[2] = z;
}

inline double
Vector::length() const
{
   return std::sqrt(_r[0]*_r[0] + _r[1]*_r[1] + _r[2]*_r[2]);
}

inline void
Vector::normalize()
{
   double len = length();
   if ( len > 0 )
      *this/=len; 
}
   
inline Vector&
Vector::operator+=(const Vector& a)
{
   this->_r[0] += a._r[0];
   this->_r[1] += a._r[1];
   this->_r[2] += a._r[2];
   return *this;
}
inline Vector&
Vector::operator-=(const Vector& a)
{
   this->_r[0] -= a._r[0];
   this->_r[1] -= a._r[1];
   this->_r[2] -= a._r[2];
   return *this;
}
inline Vector&
Vector::operator*=(const double& a)
{
   this->_r[0] *= a;
   this->_r[1] *= a;
   this->_r[2] *= a;
   return *this;
}
inline Vector&
Vector::operator/=(const double& a)
{
   this->_r[0] /= a;
   this->_r[1] /= a;
   this->_r[2] /= a;
   return *this;
}

inline Vector
Vector::operator-() const
{
   return (*this) * -1.0;
}


inline bool
Vector::operator<(const Vector& b) const
{
   return
      _r[0]<b._r[0] ||
      (_r[0]==b._r[0] &&
       (_r[1]<b._r[1] ||
        ( _r[1] == b._r[1] && _r[2] < b._r[2] )));
}

inline Vector
operator+(const Vector& a, const Vector& b)
{
   Vector c(a);
   return c+=b;
}
inline Vector
operator-(const Vector& a, const Vector& b)
{
   Vector c(a);
   return c-=b;
}
inline Vector
operator*(const Vector& a, const double& b)
{
   Vector c(a);
   return c*=b;
}
inline Vector
operator*(const double& a, const Vector& b)
{
   Vector c(b);
   return c*=a;
}
inline Vector
operator/(const Vector& a, const double& b)
{
   Vector c(a);
   return c/=b;
}
inline Vector
operator/(const double& a, const Vector& b)
{
   Vector c(b);
   return c/=a;
}
inline double
dot(const Vector& a, const Vector& b)
{
   return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline double
length(const Vector& a)
{
   return a.length();
}

inline double diffSq(const Vector& a, const Vector& b)
{
   return dot(a-b, a-b);
}

#endif // ifdef VECTOR_HH
