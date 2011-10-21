#ifndef FD_LAPLACIAN_HH
#define FD_LAPLACIAN_HH

#include "conductivity.hh"

inline double
boundaryFDLaplacianSaleheen98SumPhi(double*** &phi, Diffusion*** &diffConst, int &x, int &y, int &z)
{
   
   double SumAphi = (diffConst[x][y][z].A1  * (phi[x+1][y][z]))
                  + (diffConst[x][y][z].A2  * (phi[x][y+1][z]))
                  + (diffConst[x][y][z].A3  * (phi[x-1][y][z]))
                  + (diffConst[x][y][z].A4  * (phi[x][y-1][z]))
                  + (diffConst[x][y][z].A5  * (phi[x+1][y+1][z]))
                  + (diffConst[x][y][z].A6  * (phi[x-1][y+1][z]))
                  + (diffConst[x][y][z].A7  * (phi[x-1][y-1][z]))
                  + (diffConst[x][y][z].A8  * (phi[x+1][y-1][z]))
                  + (diffConst[x][y][z].A9  * (phi[x][y][z+1]))
                  + (diffConst[x][y][z].A10 * (phi[x][y][z-1]))
                  + (diffConst[x][y][z].A11 * (phi[x][y+1][z+1]))
                  + (diffConst[x][y][z].A12 * (phi[x][y+1][z-1]))
                  + (diffConst[x][y][z].A13 * (phi[x][y-1][z-1]))
                  + (diffConst[x][y][z].A14 * (phi[x][y-1][z+1]))
                  + (diffConst[x][y][z].A15 * (phi[x+1][y][z+1]))
                  + (diffConst[x][y][z].A16 * (phi[x-1][y][z+1]))
                  + (diffConst[x][y][z].A17 * (phi[x-1][y][z-1]))
                  + (diffConst[x][y][z].A18 * (phi[x+1][y][z-1]));
  
  double result = SumAphi - (diffConst[x][y][z].sumA * (phi[x][y][z]));
  return(result);
}


#endif
