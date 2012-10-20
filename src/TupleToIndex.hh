#ifndef TUPLE_TO_INDEX_HH
#define TUPLE_TO_INDEX_HH

#include "Long64.hh"

class TupleToIndex
{
 public:
   TupleToIndex(int nx, int ny, int nz);
   
   Long64 const operator()(int ix, int iy, int iz);
   
 private:
    Long64 nx_;  // needs to be Long64 to force 64-bit math below
    Long64 ny_;
    Long64 nz_;
};

inline TupleToIndex::TupleToIndex(int nx, int ny, int nz)
: nx_(nx), ny_(ny), nz_(nz)
{}

inline Long64
const TupleToIndex::operator()(int ix, int iy, int iz)
{
   return ix + nx_*(iy + ny_*(iz));
}

#endif
