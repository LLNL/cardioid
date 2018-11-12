#ifndef INDEX_TO_TUPLE_HH
#define INDEX_TO_TUPLE_HH

#include "Tuple.hh"
#include "Long64.hh"

class IndexToTuple
{
 public:
   IndexToTuple(int nx, int ny, int nz)
   : nx_(nx), ny_(ny), nz_(nz)
   {};

   Tuple operator()(Long64 index) const 
   {
      int x = index % nx_;
      index /= nx_;
      int y = index % ny_;
      int z = index / ny_;

      return Tuple(x, y, z);
   }
   
 private:
   int nx_;
   int ny_;
   int nz_;
};

#endif
