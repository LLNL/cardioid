#ifndef GRID3DSTENCIL_H
#define GRID3DSTENCIL_H

#include <vector>
using std::vector;

class Grid3DStencil
{
  private:

  int nstencil_;          // stencil size
  int gid_;
  int nx_, ny_, nz_;      // grid dimensions
  vector<int> nbr_gids_;  // gids of stencil points
  
  public:
  
  Grid3DStencil(int gid, int nx, int ny, int nz);
  ~Grid3DStencil(void);
  int nstencil(void) { return nstencil_; }
  int nbr_gid(int i) { return nbr_gids_[i]; }
  vector<int> nbr_gid(void) { return nbr_gids_; }
};
#endif
