#ifndef GRIDPOINT_H
#define GRIDPOINT_H
class GridPoint {

  public:

  int x,y,z;

  GridPoint(int gid, int nx, int ny, int nz)
  {
    x = gid % nx;
    gid /= nx;
    y = gid % ny;
    z = gid / ny;
  }
};
#endif
