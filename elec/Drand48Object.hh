#ifndef DRAND48_OBJECT_HH
#define DRAND48_OBJECT_HH

#include <cstdlib>

class Drand48Object
{
 public:
   Drand48Object(long seed)
   {
      srand48(seed);
      for (unsigned ii=0; ii<2000; ++ii)
         drand48();
   }
   
   int operator()(int maxval)
   {
      return int(drand48() * maxval);
   }
   
};

#endif
