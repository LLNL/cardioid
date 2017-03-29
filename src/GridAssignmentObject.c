#include "GridAssignmentObject.h"

#include <assert.h>
#include <math.h>

#include "three_algebra.h"
#include "ddcMalloc.h"

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))

static THREE_INT whichCellTuple(const GRID_ASSIGNMENT_OBJECT* this, const THREE_VECTOR r);
static int whichCell(const GRID_ASSIGNMENT_OBJECT* this, const THREE_VECTOR r);
static int tupleToIndex(const GRID_ASSIGNMENT_OBJECT* this, THREE_INT tuple);
static THREE_INT indexToTuple(const GRID_ASSIGNMENT_OBJECT* this, int index);
static double minDist2(const GRID_ASSIGNMENT_OBJECT* this, const THREE_VECTOR r, int iCell);
static void addTupleToQueue(GRID_ASSIGNMENT_OBJECT*this, THREE_INT iTuple);
static void addNbrsToQueue(GRID_ASSIGNMENT_OBJECT* this, int iCell);

/** The present implementation of GRID_ASSIGNMENT_OBJECT is judged to be
 *  sufficiently fast to meet the needs of initial assignment of 
 *  particles to domains.  The best way to speed up the code would be to
 *  more strictly limit the number of cells that are flooded by 
 *  implementing an improved distance calculation in minDist2.
 *
 *  The next best optimization possibility probably involves reducing
 *  the number of indexToTuple and tupleToIndex calculations (probably
 *  at the expense of a higher memory footprint.
*/


GRID_ASSIGNMENT_OBJECT* gao_init(int nCenters,
				 const void* centerP,
				 int stride)
{
   GRID_ASSIGNMENT_OBJECT* this = ddcMalloc(sizeof(GRID_ASSIGNMENT_OBJECT));
   this->nCenters = nCenters;
   this->centerP = centerP;
   this->stride = stride;
   this->nextCenter = ddcMalloc(nCenters*sizeof(int));
   this->floodQueue = intQueue_init(200);
   this->burnList = intQueue_init(200);
   
   // This sets the length scale of the grid cells.  The value 5 is
   // pretty arbitrary.  It could just as easily be 1 or 10.  If
   // necessary it could be made a parameter that is wired out to the
   // input deck.
   int centersPerCell = 5;

   THREE_VECTOR minCoord = *((THREE_VECTOR*) centerP);
   THREE_VECTOR maxCoord = *((THREE_VECTOR*) centerP);
   for (int ii=1; ii<nCenters; ++ii)
   {
      THREE_VECTOR iCenter = *(THREE_VECTOR*)(((char*)centerP) + ii*stride);
      minCoord.x = MIN(minCoord.x, iCenter.x);
      minCoord.y = MIN(minCoord.y, iCenter.y);
      minCoord.z = MIN(minCoord.z, iCenter.z);
      maxCoord.x = MAX(maxCoord.x, iCenter.x);
      maxCoord.y = MAX(maxCoord.y, iCenter.y);
      maxCoord.z = MAX(maxCoord.z, iCenter.z);
   }
   this->corner = minCoord;

   // It is possible that all of the centers lie on the x-, y-, or
   // z-plane.  If so, arbitrarily set the length in that direction to
   // 1.
   double lx = MAX(1, (maxCoord.x - minCoord.x));
   double ly = MAX(1, (maxCoord.y - minCoord.y));
   double lz = MAX(1, (maxCoord.z - minCoord.z));

   double x = nCenters/centersPerCell/(lx*ly*lz);
   x = pow(x, 1.0/3.0);
   this->nx = MAX(1, floor(x*lx));
   this->ny = MAX(1, floor(x*ly));
   this->nz = MAX(1, floor(x*lz));
   this->dx = lx/this->nx;
   this->dy = ly/this->ny;
   this->dz = lz/this->nz;

   int nCells = this->nx * this->ny * this->nz;
   this->grid = ddcMalloc(nCells*sizeof(GAO_CELL));
   for (int ii=0; ii<nCells; ++ii)
   {
      this->grid[ii].burned = 0;
      this->grid[ii].first = -1;
   }
   
   for (int ii=0; ii<nCenters; ++ii)
   {
      THREE_VECTOR iCenter = *(THREE_VECTOR*)(((char*)centerP) + ii*stride);
      int iCell = whichCell(this, iCenter);
      this->nextCenter[ii] = this->grid[iCell].first;
      this->grid[iCell].first = ii;
   }
   return this;
   
}

void gao_destroy(GRID_ASSIGNMENT_OBJECT* this)
{
   ddcFree(this->grid);
   intQueue_destroy(this->burnList);
   intQueue_destroy(this->floodQueue);
   ddcFree(this->nextCenter);
}

int gao_nearestCenter(GRID_ASSIGNMENT_OBJECT* this, const THREE_VECTOR r)
{
   double r2Min = 1e300;
   int minCenter = -1;
   THREE_INT iTuple = whichCellTuple(this, r);
   addTupleToQueue(this, iTuple);

   while (intQueue_size(this->floodQueue) > 0)
   {
      // pop the next cell to check
      int iCell = intQueue_pop(this->floodQueue);
      // if cell is too far away to bother continue.
      if (minDist2(this, r, iCell) > r2Min)
	 continue;
      // check all centers in this cell
      int iCenter = this->grid[iCell].first;
      while (iCenter >= 0)
      {
	 THREE_VECTOR rCenter =
	    *(THREE_VECTOR*)(((char*)this->centerP) + iCenter*this->stride);
	 double r2 = DIFFSQ(r, rCenter);
         if (r2 == r2Min)
            minCenter = MIN(minCenter, iCenter);
	 if (r2 < r2Min)
	 {
	    r2Min = r2;
	    minCenter = iCenter;
	 }
	 iCenter = this->nextCenter[iCenter];
      }
      // push any unburned nbrs to queue, burnList.  Mark as burned
      addNbrsToQueue(this, iCell);
   }

   while (intQueue_size(this->burnList) > 0)
      this->grid[intQueue_pop(this->burnList)].burned = 0;
   
   assert(minCenter >= 0);
   return minCenter;
}

      
THREE_INT whichCellTuple(const GRID_ASSIGNMENT_OBJECT* this, const THREE_VECTOR r)
{
   THREE_INT iTuple;
   iTuple.x = (r.x-this->corner.x)/this->dx;
   iTuple.y = (r.y-this->corner.y)/this->dy;
   iTuple.z = (r.z-this->corner.z)/this->dz;
   iTuple.x = MAX(0, iTuple.x);
   iTuple.y = MAX(0, iTuple.y);
   iTuple.z = MAX(0, iTuple.z);  
   iTuple.x = MIN(this->nx-1, iTuple.x);
   iTuple.y = MIN(this->ny-1, iTuple.y);
   iTuple.z = MIN(this->nz-1, iTuple.z);

   return iTuple;
}

int whichCell(const GRID_ASSIGNMENT_OBJECT* this, const THREE_VECTOR r)
{
   return tupleToIndex(this, whichCellTuple(this, r));
}

int tupleToIndex(const GRID_ASSIGNMENT_OBJECT* this, THREE_INT tuple)
{
   return tuple.x + this->nx * (tuple.y + this->ny*tuple.z);
}

THREE_INT indexToTuple(const GRID_ASSIGNMENT_OBJECT* this, int index)
{
   THREE_INT iTuple;
   iTuple.x = index % this->nx;
   index /= this->nx;
   iTuple.y = index % this->ny;
   iTuple.z = index / this->ny;
   return iTuple;
}

/** Finds a lower bound of the squared distance from the point r to the
 * cell with index iCell.  As presently implemented this calculation is
 * very conservative.  We could set a larger lower bound by considering
 * the location of the particle within the cell in which it lies.  */
double minDist2(const GRID_ASSIGNMENT_OBJECT* this, const THREE_VECTOR r, int iCell)
{
   THREE_INT ir = whichCellTuple(this, r);
   THREE_INT iTuple = indexToTuple(this, iCell);
   
   double rx = this->dx*(abs(iTuple.x - ir.x) - 1); rx = MAX(0, rx);
   double ry = this->dy*(abs(iTuple.y - ir.y) - 1); ry = MAX(0, ry);
   double rz = this->dz*(abs(iTuple.z - ir.z) - 1); rz = MAX(0, rz);

   return rx*rx + ry*ry + rz*rz;
}

void addTupleToQueue(GRID_ASSIGNMENT_OBJECT*this, THREE_INT iTuple)
{
   int index = tupleToIndex(this, iTuple);
   if (this->grid[index].burned != 0)
      return;
   intQueue_push(this->floodQueue, index);
   intQueue_push(this->burnList, index);
   this->grid[index].burned = 1;
}

void addNbrsToQueue(GRID_ASSIGNMENT_OBJECT* this, int iCell)
{
   THREE_INT iTuple = indexToTuple(this, iCell);
   iTuple.x += 1; if (iTuple.x < this->nx) addTupleToQueue(this, iTuple);
   iTuple.x -= 2; if (iTuple.x >= 0)       addTupleToQueue(this, iTuple);
   iTuple.x += 1;
   
   iTuple.y += 1; if (iTuple.y < this->ny) addTupleToQueue(this, iTuple);
   iTuple.y -= 2; if (iTuple.y >= 0)       addTupleToQueue(this, iTuple);
   iTuple.y += 1;

   iTuple.z += 1; if (iTuple.z < this->nz) addTupleToQueue(this, iTuple);
   iTuple.z -= 2; if (iTuple.z >= 0)       addTupleToQueue(this, iTuple);
   iTuple.z += 1;
}

