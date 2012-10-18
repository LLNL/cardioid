#ifndef GRID3DSTENCIL_H
#define GRID3DSTENCIL_H

#include <vector>
#include "Long64.hh"
using std::vector;

class Grid3DStencil
{
 public:
   
   Grid3DStencil(Long64 gid, int nx, int ny, int nz);
   
   int nStencil(void) const { return nbrGids_.size(); }
   Long64 nbrGid(int i) const { return nbrGids_[i]; }
   const Long64& operator[](int i) const { return nbrGids_[i]; }
   const vector<Long64>& nbrGid(void) const { return nbrGids_; }

 private:
   
   vector<Long64> nbrGids_;   // gids of stencil points
   
};
#endif
