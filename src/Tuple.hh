#ifndef TUPLE_HH
#define TUPLE_HH

class Tuple
{
 public:
   Tuple(int ix, int iy, int iz);

   const int& x() const;
   const int& y() const;
   const int& z() const;

   Tuple& operator-=(const Tuple& a);
   Tuple& operator+=(const Tuple& a);
   
 private:
   int ix_;
   int iy_;
   int iz_;
};

Tuple operator-(const Tuple& a, const Tuple& b);
Tuple operator+(const Tuple& a, const Tuple& b);




inline Tuple::Tuple(int ix, int iy, int iz)
: ix_(ix), iy_(iy), iz_(iz)
{
}

inline const int& Tuple::x() const
{
   return ix_;
}

inline const int& Tuple::y() const 
{
   return iy_;
}

inline const int& Tuple::z() const
{
   return iz_;
}

inline Tuple& Tuple::operator-=(const Tuple& a)
{
   ix_ -= a.ix_;
   iy_ -= a.iy_;
   iz_ -= a.iz_;
}

inline Tuple& Tuple::operator+=(const Tuple& a)
{
   ix_ += a.ix_;
   iy_ += a.iy_;
   iz_ += a.iz_;
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

#endif
