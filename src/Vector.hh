#ifndef VECTOR_HH
#define VECTOR_HH

#include <iosfwd>

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

inline double diffSq(const Vector& a, const Vector& b)
{
   return dot(a-b, a-b);
}



#endif // ifdef VECTOR_HH
