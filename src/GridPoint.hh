#ifndef GRIDPOINT_H
#define GRIDPOINT_H
#include "Long64.hh"
class GridPoint {

  public:

  int x,y,z;

  GridPoint(Long64 gid, int nx, int ny, int nz)
  {
    x = gid % nx;
    gid /= nx;
    y = gid % ny;
    z = gid / ny;
  }
};
#endif
