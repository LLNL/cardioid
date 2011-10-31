#ifndef INDEX_TO_THREE_VECTOR_HH
#define INDEX_TO_THREE_VECTOR_HH

#include "three_algebra.h"
#include "Long64.hh"

class IndexToThreeVector
{
 public:
   IndexToThreeVector(int nx, int ny, int nz)
   : nx_(nx), ny_(ny), nz_(nz)
   {};

   THREE_VECTOR operator()(Long64 index) const
   {
      THREE_VECTOR r;
      r.x = index % nx_;
      index /= nx_;
      r.y = index % ny_;
      r.z = index / ny_;

      return r;
   }
   
 private:
   int nx_;
   int ny_;
   int nz_;
};

#endif
