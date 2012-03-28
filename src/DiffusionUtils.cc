#include "DiffusionUtils.hh"
#include "Anatomy.hh"
#include "LocalGrid.hh"
#include "Tuple.hh"

#include <algorithm>

using namespace std;

/** We want to find the boundingBox such that any stencil point of any
 *  local atom is in the box.  It is not sufficient merely to iterate all
 *  of the local and remote atoms and find the maximum extent.  There
 *  may be local cells that are on the outer or inner walls of the
 *  heart.  Such cells will have no remote cells to satisfy their
 *  stencil.  Therefore, the safe bet is to iterate the local cells and
 *  add the stencil size in each direction.
 */
LocalGrid DiffusionUtils::findBoundingBox(const Anatomy& anatomy)
{
   assert(anatomy.nLocal() > 0);
   Tuple globalTuple = anatomy.globalTuple(0);
   int xMin = globalTuple.x();
   int yMin = globalTuple.y();
   int zMin = globalTuple.z();
   int xMax = globalTuple.x();
   int yMax = globalTuple.y();
   int zMax = globalTuple.z();
   
   for (unsigned ii=1; ii<anatomy.nLocal(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      xMin = min(xMin, globalTuple.x());
      yMin = min(yMin, globalTuple.y());
      zMin = min(zMin, globalTuple.z());
      xMax = max(xMax, globalTuple.x());
      yMax = max(yMax, globalTuple.y());
      zMax = max(zMax, globalTuple.z());
   }
   
   int stencilSize = 1;
   
   int nx = 2*stencilSize + xMax - xMin + 1;
   int ny = 2*stencilSize + yMax - yMin + 1;
   int nz = 2*stencilSize + zMax - zMin + 1;
   xMin -= stencilSize;
   yMin -= stencilSize;
   zMin -= stencilSize;

   return LocalGrid(nx, ny, nz, xMin, yMin, zMin);
};

LocalGrid DiffusionUtils::findBoundingBox_simd(const Anatomy& anatomy)
{
   assert(anatomy.nLocal() > 0);
   Tuple globalTuple = anatomy.globalTuple(0);
   int xMin = globalTuple.x();
   int yMin = globalTuple.y();
   int zMin = globalTuple.z();
   int xMax = globalTuple.x();
   int yMax = globalTuple.y();
   int zMax = globalTuple.z();
   
   for (unsigned ii=1; ii<anatomy.nLocal(); ++ii)
   {
      Tuple globalTuple = anatomy.globalTuple(ii);
      xMin = min(xMin, globalTuple.x());
      yMin = min(yMin, globalTuple.y());
      zMin = min(zMin, globalTuple.z());
      xMax = max(xMax, globalTuple.x());
      yMax = max(yMax, globalTuple.y());
      zMax = max(zMax, globalTuple.z());
   }
   
   int stencilSize = 1;
   
   int nx = 2*stencilSize + xMax - xMin + 1;
   int ny = 2*stencilSize + yMax - yMin + 1;
   int nz = 2*stencilSize + zMax - zMin + 1;
   xMin -= stencilSize;
   yMin -= stencilSize;
   zMin -= stencilSize;

   nz += (nz%4==0) ? 0 : 4 - (nz%4);
   
   return LocalGrid(nx, ny, nz, xMin, yMin, zMin);
}
