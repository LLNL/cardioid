#ifndef INDEX_TO_VECTOR_HH
#define INDEX_TO_VECTOR_HH

#include "Vector.hh"
#include "Long64.hh"

class IndexToVector
{
 public:
   IndexToVector(int nx, int ny, int nz)
   : nx_(nx), ny_(ny), nz_(nz)
   {};

   Vector operator()(Long64 index) const
   {
      Vector r;
      r[0] = index % nx_;
      index /= nx_;
      r[1] = index % ny_;
      r[2] = index / ny_;

      return r;
   }
   
 private:
   int nx_;
   int ny_;
   int nz_;
};

#endif
