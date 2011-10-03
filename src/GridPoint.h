#ifndef GRIDPOINT_H
#define GRIDPOINT_H
class GridPoint {

  public:

  int x,y,z;

  GridPoint(int gid, int nx, int ny, int nz)
  {
    int nxy = nx*ny;
    z = gid/nxy;
    int xy0 = gid;
    while (xy0 >= nxy) xy0 -= nxy;
    y = xy0/nx;
    x = xy0 - y*nx;
  }
};
#endif
