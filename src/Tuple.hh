#ifndef TUPLE_HH
#define TUPLE_HH

class Tuple
{
 public:
   Tuple(int ix, int iy, int iz) : ix_(ix), iy_(iy), iz_(iz){}

   const int& x() const {return ix_;}
   const int& y() const {return iy_;}
   const int& z() const {return iz_;}

   int& x() {return ix_;}
   int& y() {return iy_;}
   int& z() {return iz_;}

   Tuple& operator-=(const Tuple& a);
   Tuple& operator+=(const Tuple& a);
   
 private:
   int ix_;
   int iy_;
   int iz_;
};

Tuple operator-(const Tuple& a, const Tuple& b);
Tuple operator+(const Tuple& a, const Tuple& b);
bool operator==(const Tuple& a, const Tuple& b);



inline Tuple& Tuple::operator-=(const Tuple& a)
{
   ix_ -= a.ix_;
   iy_ -= a.iy_;
   iz_ -= a.iz_;
   return *this;
}

inline Tuple& Tuple::operator+=(const Tuple& a)
{
   ix_ += a.ix_;
   iy_ += a.iy_;
   iz_ += a.iz_;
   return *this;
}

inline Tuple operator-(const Tuple& a, const Tuple& b)
{
   Tuple c(a);
   return c-=b;
}

inline Tuple operator+(const Tuple& a, const Tuple& b)
{
   Tuple c(a);
   return c+=b;
}

inline bool operator==(const Tuple& a, const Tuple& b)
{
   return (a.x() == b.x() && a.y() == b.y() && a.z() == b.z());
}

#endif
