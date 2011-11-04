#include "Vector.hh"

#include <iostream>

using std::ostream;

ostream&
operator<<(ostream& os, const Vector& r)
{
   os << "(" << r[0] << ", " << r[1] << ", " << r[2] << ")";
   return os;
}
